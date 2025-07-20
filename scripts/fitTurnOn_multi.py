import argparse
import os
import sys
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import copy
from scipy import interpolate
from scipy.ndimage import gaussian_filter1d

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

sys.path.insert(0, 'Common/python')
from RootObjects import Histogram, Graph

from array import array

# --------------
# Example Command:
# pip install sklearn on lxplus
# python3 scripts/fitTurnOn_multi.py --input TurnOnDeepTau/TurnOnDeepTau.root --output fitTurnOnDeepTau --decay_modes DeepTau
# python3 scripts/fitTurnOn_multi.py --input TurnOnPNet/TurnOnPNet.root --output fitTurnOnPNet --decay_modes PNet
# --------------


parser = argparse.ArgumentParser(description='Fit turn-on curves.')
parser.add_argument('--input', required=True, type=str, help="ROOT file with turn-on curves")
parser.add_argument('--output', required=True, type=str, help="output file prefix")
parser.add_argument('--channels', required=False, type=str, default='etau,mutau,ditau,ditaujet', help="channels to process")
parser.add_argument('--decay-modes', required=False, type=str, default='all,0,1,2,10,11', help="decay modes to process")
parser.add_argument('--decay_modes', required=True, type=str, default='DeepTau', choices=['DeepTau', 'PNet'], help="Type of decay modes to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight',
                    help="working points to process")
args = parser.parse_args()

def MinTarget(dy, eff):
    y = np.cumsum(dy)
    return np.sum(((eff.y - y) / (eff.y_error_high + eff.y_error_low)) ** 2)

class FitResults:
    def __init__(self, eff, x_pred):
        kernel_high = ConstantKernel()
        kernel_low = ConstantKernel() * Matern(nu=1, length_scale_bounds=(10, 100), length_scale=20)
        N = eff.x.shape[0]
        dy_init = np.gradient(eff.y, eff.x)
        dy_init = np.clip(dy_init, 0, 1)
        res = scipy.optimize.minimize(MinTarget, dy_init, args=(eff,), bounds = [ [0, 1] ] * N,
                                      options={"maxfun": int(1e6)})
        if not res.success:
            print(res)
            raise RuntimeError("Unable to prefit")

        eff = copy.deepcopy(eff)
        new_y = np.cumsum(res.x)
        delta = eff.y - new_y
        eff.y_error_low = np.sqrt(eff.y_error_low ** 2 + delta ** 2)
        eff.y_error_high = np.sqrt(eff.y_error_high ** 2 + delta ** 2)
        eff.y = new_y
        yerr = np.maximum(eff.y_error_low, eff.y_error_high)

        # self.pt_start_flat = eff.x[-1]
        # best_chi2_ndof = math.inf
        # for n in range(1, N):
        #     flat_eff, residuals, _, _, _ = np.polyfit(eff.x[N-n-1:], eff.y[N-n-1:], 0, w=1/yerr[N-n-1:], full=True)
        #     chi2_ndof = residuals[0] / n
        #     #print(n, chi2_ndof)
        #     if (chi2_ndof > 0 and chi2_ndof < best_chi2_ndof) or eff.x[N-n-1] + eff.x_error_high[N-n-1] >= 100:
        #         self.pt_start_flat = eff.x[N-n-1]
        #         best_chi2_ndof = chi2_ndof
        # if best_chi2_ndof > 20:
        #     print("Unable to determine the high pt region")
        #     self.pt_start_flat = eff.x[-1]

        self.pt_start_flat = eff.x[-1]
        best_chi2_ndof = math.inf
        # for n in range(2, N):  # require at least 2 points
        #     start_idx = N - n
        #     x_fit = eff.x[start_idx:]
        #     y_fit = eff.y[start_idx:]
        #     yerr_fit = yerr[start_idx:]

        #     # Force plateau to end value
        #     flat_val = eff.y[-1]
        #     residuals = np.sum(((y_fit - flat_val) / yerr_fit) ** 2)
        #     chi2_ndof = residuals / (len(y_fit) or 1)

        #     # Require that the slope in this region is small
        #     slope = np.gradient(y_fit, x_fit)
        #     if np.max(np.abs(slope)) > 1e-5:
        #         continue

        #     # Penalty to encourage later start for plateau
        #     pt_penalty = 0.01 * max(0, eff.x[-1] - eff.x[start_idx])
        #     chi2_penalized = chi2_ndof + pt_penalty

        #     # Accept if it's the best or very late plateau
        #     if chi2_penalized < best_chi2_ndof or (eff.x[start_idx] + eff.x_error_high[start_idx] >= 100):
        #         self.pt_start_flat = eff.x[start_idx]
        #         best_chi2_ndof = chi2_penalized
        for n in range(1, N):
            flat_eff, residuals, _, _, _ = np.polyfit(eff.x[N-n-1:], eff.y[N-n-1:], 0, w=1/yerr[N-n-1:], full=True)
            chi2_ndof = residuals[0] / n
            #print(n, chi2_ndof)
            if (chi2_ndof > 0 and chi2_ndof < best_chi2_ndof) or eff.x[N-n-1] + eff.x_error_high[N-n-1] >= 100:
                self.pt_start_flat = eff.x[N-n-1]
                best_chi2_ndof = chi2_ndof
        if best_chi2_ndof > 20:
            print("Unable to determine the high pt region")
            self.pt_start_flat = eff.x[-1]

        low_pt = eff.x <= self.pt_start_flat
        high_pt = eff.x >= self.pt_start_flat

        self.gp_high = GaussianProcessRegressor(kernel=kernel_high, alpha=yerr[high_pt] ** 2, n_restarts_optimizer=10)
        self.gp_high.fit(np.atleast_2d(eff.x[high_pt]).T, eff.y[high_pt])
        self.gp_low = GaussianProcessRegressor(kernel=kernel_low, alpha=np.append([0], yerr[low_pt] ** 2),
                                               n_restarts_optimizer=10)
        self.gp_low.fit(np.atleast_2d(np.append([10], eff.x[low_pt])).T, np.append([0], eff.y[low_pt]))

        self.y_pred, sigma_pred = self.Predict(x_pred)

        sigma_orig = np.zeros(N)
        for n in range(N):
            idx = np.argmin(abs(x_pred - eff.x[n]))
            sigma_orig[n] = sigma_pred[idx]

        interp_kind = 'linear'
        sp = interpolate.interp1d(eff.x, sigma_orig, kind=interp_kind, fill_value="extrapolate")
        sigma_interp = sp(x_pred)
        max_unc = 0.05 / math.sqrt(2)
        sigma_pred, = self.ApplyStep(x_pred, [ [ sigma_pred, sigma_interp ] ], eff.x[0], eff.x[-1] )
        outer_trend = np.minimum(np.ones(x_pred.shape[0]), (x_pred - eff.x[-1]) / eff.x[-1])
        outer_sigma = np.maximum(sigma_pred, sigma_pred + (max_unc - sigma_pred) * outer_trend )
        self.sigma_pred = np.where(x_pred < eff.x[-1], sigma_pred, outer_sigma )

    def Predict(self, x_pred):
        y_pred_high, sigma_high = self.gp_high.predict(np.atleast_2d(x_pred).T, return_std=True)
        y_pred_low, sigma_low = self.gp_low.predict(np.atleast_2d(x_pred).T, return_std=True)
        return self.ApplyStep(x_pred, [ [y_pred_low, y_pred_high], [sigma_low, sigma_high] ], self.pt_start_flat)

    def ApplyStep(self, x_pred, functions, x0, x1 = None):
        step = (np.tanh(0.1*(x_pred - x0)) + 1) / 2
        if x1 is not None:
            step *= (np.tanh(0.1*(x1 - x_pred)) + 1) / 2
        step = np.where(step > 0.999, 1, step)
        step = np.where(step < 0.001, 0, step)
        results = []
        for fn in functions:
            results.append(fn[0] * (1-step) + fn[1] * step)
        return tuple(results)

channels = args.channels.split(',')
if args.decay_modes == 'PNet':
    decay_modes = [ 'all', '0', '1', '2', '10', '11', '1011']
elif args.decay_modes == 'DeepTau':
    decay_modes = [ 'all', '0', '1', '10', '11', '1011']
working_points = args.working_points.split(',')
ch_validity_thrs = { 'etau': 35, 'mutau': 32, 'ditau': 40, 'ditaujet': 40, }

file = ROOT.TFile(args.input, 'READ')
output_dir = os.path.join(os.getcwd(), args.output)
os.makedirs(output_dir, exist_ok=True)

output_file_path = os.path.join(output_dir, f'fitTurnOn{args.decay_modes}')
print('Output file will be saved to {}'.format(output_file_path))
output_file = ROOT.TFile('{}.root'.format(output_file_path), 'RECREATE', '', ROOT.RCompressionSetting.EDefaults.kUseSmallest)

for channel in channels:
    with PdfPages('{}_{}.pdf'.format(output_file_path, channel)) as pdf:
        for wp in working_points:
            for dm in decay_modes:
                print('Processing {} {} WP DM = {}'.format(channel, wp, dm))
                dm_label = '_dm{}'.format(dm) if dm != 'all' else ''
                name_pattern = '{{}}_{}_{}{}_fit_eff'.format(channel, wp, dm_label)
                dm_label = '_dm'+ dm if len(dm) > 0 else ''
                eff_data_root = file.Get(name_pattern.format('data'))
                eff_mc_root = file.Get(name_pattern.format('mc'))
                eff_data_orig = Graph(root_graph=eff_data_root)
                eff_mc_orig = Graph(root_graph=eff_mc_root)
                pred_step = 0.1
                #x_low = min(eff_data.x[0] - eff_data.x_error_low[0], eff_mc.x[0] - eff_mc.x_error_low[0])
                #x_high = max(eff_data.x[-1] + eff_data.x_error_high[-1], eff_mc.x[-1] + eff_mc.x_error_high[-1])
                x_low, x_high = 20, 1000
                x_pred = np.arange(x_low, x_high + pred_step / 2, pred_step)


                def get_common_bins(graph_data, graph_mc, rel_unc_threshold=0.3, min_width=2.0, max_merge=10):
                    """
                    Merge bins from high to low pT. Ensures no empty bins in either data or MC.
                    Returns final merged bin edges (to use identically for both).
                    """
                    def extract_bins(graph):
                        bins = []
                        for i in range(graph.GetN()):
                            x = graph.GetX()[i]
                            ex_low = graph.GetErrorXlow(i)
                            ex_high = graph.GetErrorXhigh(i)
                            y = graph.GetY()[i]
                            ey = 0.5 * (graph.GetErrorYlow(i) + graph.GetErrorYhigh(i))
                            rel_unc = ey / y if y > 0 else float('inf')
                            bins.append({'low': x - ex_low, 'high': x + ex_high, 'rel_unc': rel_unc})
                        return bins

                    bins_data = extract_bins(graph_data)
                    bins_mc = extract_bins(graph_mc)
                    combined = sorted(bins_data + bins_mc, key=lambda b: b['low'])

                    # Full sorted unique edges from both
                    edges = sorted(set([b['low'] for b in combined] + [b['high'] for b in combined]))
                    
                    merged_edges = [edges[-1]]
                    i = len(edges) - 2
                    while i >= 0:
                        bin_low = edges[i]
                        bin_high = merged_edges[0]
                        merge_count = 0

                        while True:
                            # Check if bin [bin_low, bin_high] is empty in either graph
                            empty_data = True
                            for j in range(graph_data.GetN()):
                                x = graph_data.GetX()[j]
                                if bin_low <= x < bin_high:
                                    empty_data = False
                                    break
                            empty_mc = True
                            for j in range(graph_mc.GetN()):
                                x = graph_mc.GetX()[j]
                                if bin_low <= x < bin_high:
                                    empty_mc = False
                                    break

                            # Check combined relative uncertainty
                            in_data = [b['rel_unc'] for b in bins_data if b['low'] >= bin_low and b['high'] <= bin_high]
                            in_mc = [b['rel_unc'] for b in bins_mc if b['low'] >= bin_low and b['high'] <= bin_high]
                            all_uncs = in_data + in_mc
                            mean_unc = np.mean(all_uncs) if all_uncs else float('inf')
                            width = bin_high - bin_low

                            # Stop merging if: has data & mc & good uncertainty & width, or max_merge reached
                            if (not empty_data and not empty_mc and mean_unc < rel_unc_threshold and width >= min_width) or \
                            merge_count >= max_merge or i == 0:
                                break

                            # Merge backward
                            merge_count += 1
                            i -= 1
                            bin_low = edges[max(0, i)]

                        merged_edges.insert(0, bin_low)
                        i -= 1

                    return np.array(merged_edges)


                # Interpolate both onto common bin centers if needed
                bin_edges = get_common_bins(eff_data_root, eff_mc_root)
                bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

                def rebin_graph(original_graph, bin_edges):
                    """
                    Rebins a Graph into fixed bins defined by bin_edges.
                    No merging logic — assumes bin_edges are already merged.
                    """
                    x_vals = []
                    y_vals = []
                    x_err_low = []
                    x_err_high = []
                    y_err_low = []
                    y_err_high = []

                    for i in range(len(bin_edges) - 1):
                        x_low = bin_edges[i]
                        x_high = bin_edges[i+1]

                        in_bin = (original_graph.x >= x_low) & (original_graph.x < x_high)
                        if not np.any(in_bin):
                            continue

                        x_bin = original_graph.x[in_bin]
                        y_bin = original_graph.y[in_bin]
                        y_err_low_bin = original_graph.y_error_low[in_bin]
                        y_err_high_bin = original_graph.y_error_high[in_bin]

                        weights = 1. / np.maximum(1e-12, (0.5 * (y_err_low_bin + y_err_high_bin)) ** 2)
                        weighted_avg = np.average(y_bin, weights=weights)
                        avg_err = np.sqrt(1. / np.sum(weights))

                        center = 0.5 * (x_low + x_high)
                        width = 0.5 * (x_high - x_low)

                        x_vals.append(center)
                        y_vals.append(weighted_avg)
                        x_err_low.append(width)
                        x_err_high.append(width)
                        y_err_low.append(avg_err)
                        y_err_high.append(avg_err)

                    graph_rebinned = ROOT.TGraphAsymmErrors(
                        len(x_vals),
                        array('d', x_vals),
                        array('d', y_vals),
                        array('d', x_err_low),
                        array('d', x_err_high),
                        array('d', y_err_low),
                        array('d', y_err_high),
                    )

                    return Graph(root_graph=graph_rebinned)

                eff_data = rebin_graph(eff_data_orig, bin_edges)
                eff_mc = rebin_graph(eff_mc_orig, bin_edges)

                eff_data_fitted = FitResults(eff_data, x_pred)
                eff_mc_fitted = FitResults(eff_mc, x_pred)

                sf = eff_data_fitted.y_pred / eff_mc_fitted.y_pred
                sf_sigma = np.sqrt( (eff_data_fitted.sigma_pred / eff_mc_fitted.y_pred) ** 2 \
                         + (eff_data_fitted.y_pred / (eff_mc_fitted.y_pred ** 2) * eff_mc_fitted.sigma_pred ) ** 2 )

                fig, (ax, ax_ratio) = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                                                           gridspec_kw = {'height_ratios':[2, 1]})
                mc_color = 'g'
                data_color = 'k'
                trans = 0.3

                # test by botao
                # print("mc low previous: {}".format(eff_mc.x_error_low))
                # print("data low previous: {}".format(eff_data.x_error_low))
                # count = 0
                # for _i,_obj in enumerate(eff_data.x_error_low) :
                #     if _obj < 0:
                #         _obj == abs(_obj)
                # for _i,_obj in enumerate(eff_data.x_error_high) :
                #     if _obj < 0:
                #         _obj == - _obj
                # count = 0
                # for _i,_obj in enumerate(eff_mc.x_error_low) :
                #     if _obj < 0:
                #         print("BUGBUGBUGBUG!")
                #         count = 1
                #         _obj == abs(_obj)
                # if count == 1:
                #     continue
                # for _i,_obj in enumerate(eff_mc.x_error_high) :
                #     if _obj < 0:
                #         _obj == - _obj
                # print("mc low: {}".format(eff_mc.x_error_low))
                # print("mc high: {}".format(eff_mc.x_error_high))
                # print("data low: {}".format(eff_data.x_error_low))
                # print("data high: {}".format(eff_data.x_error_high))
                # end test
                plt_data = ax.errorbar(eff_data.x, eff_data.y, xerr=(abs(eff_data.x_error_low), abs(eff_data.x_error_high)),
                                       yerr=(eff_data.y_error_low, eff_data.y_error_high), fmt=data_color+'.',
                                       markersize=5)
                plt_mc = ax.errorbar(eff_mc.x, eff_mc.y, xerr=(abs(eff_mc.x_error_low), abs(eff_mc.x_error_high)),
                                     yerr=(eff_mc.y_error_low, eff_mc.y_error_high), fmt=mc_color+'.', markersize=5)

                plt_data_fitted = ax.plot(x_pred, eff_data_fitted.y_pred, data_color+'--')
                ax.fill(np.concatenate([x_pred, x_pred[::-1]]),
                        np.concatenate([eff_data_fitted.y_pred - eff_data_fitted.sigma_pred,
                                       (eff_data_fitted.y_pred + eff_data_fitted.sigma_pred)[::-1]]),
                        alpha=trans, fc=data_color, ec='None')

                plt_mc_fitted = ax.plot(x_pred, eff_mc_fitted.y_pred, mc_color+'--')
                ax.fill(np.concatenate([x_pred, x_pred[::-1]]),
                        np.concatenate([eff_mc_fitted.y_pred - eff_mc_fitted.sigma_pred,
                                       (eff_mc_fitted.y_pred + eff_mc_fitted.sigma_pred)[::-1]]),
                        alpha=trans, fc=mc_color, ec='None')

                ax_ratio.plot(x_pred, sf, 'b--')
                ax_ratio.fill(np.concatenate([x_pred, x_pred[::-1]]),
                              np.concatenate([sf - sf_sigma, (sf + sf_sigma)[::-1]]),
                              alpha=trans, fc='b', ec='None')

                title = "Turn-ons for {} trigger with {} DeepTau VSjet".format(channel, wp)
                if dm != 'all':
                    title += " for DM={}".format(dm)
                else:
                    title += " for all DMs"
                ax.set_title(title, fontsize=16)
                ax.set_ylabel("Efficiency", fontsize=12)
                ax.set_ylim([ 0., 1.1 ])
                ax.set_xlim([ 20, min(200, plt.xlim()[1]) ])

                ax_ratio.set_xlabel("$p_T$ (GeV)", fontsize=12)
                ax_ratio.set_ylabel("Data/MC SF", fontsize=12)
                ax_ratio.set_ylim([0.5, 1.49])

                validity_plt = ax.plot( [ ch_validity_thrs[channel] ] * 2, ax.get_ylim(), 'r--' )
                ax_ratio.plot( [ ch_validity_thrs[channel] ] * 2, ax_ratio.get_ylim(), 'r--' )

                ax.legend([ plt_data, plt_mc, plt_data_fitted[0], plt_mc_fitted[0], validity_plt[0] ],
                          [ "Data", "MC", "Data fitted", "MC fitted", "Validity range"], fontsize=12, loc='lower right')


                plt.subplots_adjust(hspace=0)
                pdf.savefig(bbox_inches='tight')
                plt.close()

                out_name_pattern = '{{}}_{}_{}{}_{{}}'.format(channel, wp, dm_label)
                output_file.WriteTObject(eff_data_root, out_name_pattern.format('data', 'eff'), 'Overwrite')
                output_file.WriteTObject(eff_mc_root, out_name_pattern.format('mc', 'eff'), 'Overwrite')
                eff_data_fitted_hist = Histogram.CreateTH1(eff_data_fitted.y_pred, [x_low, x_high],
                                                           eff_data_fitted.sigma_pred, fixed_step=True)
                eff_mc_fitted_hist = Histogram.CreateTH1(eff_mc_fitted.y_pred, [x_low, x_high],
                                                         eff_mc_fitted.sigma_pred, fixed_step=True)
                sf_fitted_hist = eff_data_fitted_hist.Clone()
                sf_fitted_hist.Divide(eff_mc_fitted_hist)
                output_file.WriteTObject(eff_data_fitted_hist, out_name_pattern.format('data', 'fitted'), 'Overwrite')
                output_file.WriteTObject(eff_mc_fitted_hist, out_name_pattern.format('mc', 'fitted'), 'Overwrite')
                output_file.WriteTObject(sf_fitted_hist, out_name_pattern.format('sf', 'fitted'), 'Overwrite')

output_file.Close()
print('All done.')
os._exit(0)
