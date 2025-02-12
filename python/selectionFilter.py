# write a Module for selection filter

import sys; sys.path.append('python')

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from summaryProducer import *
from tupleProducer import *

import ROOT as R
import math

R.PyConfig.IgnoreCommandLineOptions = True

class selectionFilter(Module):
    def __init__(self, isMC, era):
        self.isMC = isMC
        self.era = era
        self.nanoVer = 11

        # cutflow hist
        self.cutflow_hist = R.TH1F('pre_selection','pre_selection',20,0,20)
        self.cutflow_assigned = 0
        pass
    
    def beginJob(self):
        pass
    
    def endJob(self):
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.cutflow_hist.Write()
        pass

    def fill_cut(self, cutnm):
        _ibin = self.cutflow_hist.GetXaxis().FindBin(cutnm)
        if _ibin == -1:
            self.cutflow_assigned += 1
            self.cutflow_hist.GetXaxis().SetBinLabel(self.cutflow_assigned, cutnm)
            self.cutflow_hist.SetBinContent( self.cutflow_assigned, 1 )
        else:
            self.cutflow_hist.SetBinContent( _ibin, self.cutflow_hist.GetBinContent(_ibin) + 1 )

    # function for calculate delta phi
    def calc_dphi(self, lepton_v4, met):    # object : lepton, met
        dphi = lepton_v4.Phi() - met.Phi()
        while dphi<-math.pi:
            dphi+=2*math.pi
        while dphi>math.pi:
            dphi-=2*math.pi
        return dphi
    
    # function for calculate MT
    def calc_MT(self, lepton_v4, met):
        delta_phi = self.calc_dphi(lepton_v4, met)
        return math.sqrt( 2.0 * lepton_v4.Pt() * met.Pt() * (1.0 - math.cos(delta_phi)) )

    def calc_met(self, met):
        _met_v = R.TLorentzVector()
        _met_v.SetPtEtaPhiM(met.pt, 0, met.phi, 0)
        return _met_v

    def find_signal_muon(self,coll_muons):
        muon_sel = []
        for _m, _muon in enumerate(coll_muons):
            _mu_v4 = _muon.p4()
            if (_muon.mediumId) and (_mu_v4.Pt() > 24) and (abs(_mu_v4.Eta()) < 2.1) :
                muon_sel.append( (_mu_v4, _muon, _m) )
        if (self.nanoVer == 10 or self.nanoVer == 11):
            muon_sel.sort(key=lambda x:(x[1].pfRelIso04_all, x[0].Pt()) ,reverse=True)
        else:
            muon_sel.sort(key=lambda x:(x[1].miniIsoId, x[0].Pt()) ,reverse=True)
        return muon_sel

    def analyze(self, event):
        # process event, return True (go to next module) or False (fail, go to next event)
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        taus = Collection(event, "Tau")
        met = Object(event, "MET")
        flag = Object(event, "Flag")
        
        signalMuon = None
        signalMuon_num = None
        find_signal_muon = False
        has_other_muon = False
        has_ele = False
        MT_Cut = False
        has_bjet = False
        pass_MET_filter = True
        mtCut = 0  #100
        btagThreshold = 0  #0.8
        
        lt_muon_sel = []
        lt_muon_sel = self.find_signal_muon(muons)
        if len(lt_muon_sel) >= 1:
            signalMuon_v4 = lt_muon_sel[0][0]
            signalMuon_num = lt_muon_sel[0][2]
            find_signal_muon = True

        # Find signal muon
        '''
        lt_muon_sel = []
        for _m, _muon in enumerate(muons):
            _mu_v4 = _muon.p4()
            if (_muon.mediumId) and (_mu_v4.Pt() > 24) and (abs(_mu_v4.Eta()) < 2.1) :
                lt_muon_sel.append( (_mu_v4, _muon, _m) )
        lt_muon_sel.sort(key=lambda x:(x[1].miniIsoId, x[0].Pt()) ,reverse=True)
        if len(lt_muon_sel) >= 1:
            signalMuon = lt_muon_sel[0][0]
            signalMuon_num = lt_muon_sel[0][2]
            print("sele muon: ",signalMuon.Pt())
            find_signal_muon = True
        '''

        # Apply third lepton veto
        for _m, _muon in enumerate(muons):
            _mu_v4 = _muon.p4()
            if (_m != signalMuon_num) and (_mu_v4.Pt() > 10) and (abs(_mu_v4.Eta()) < 2.4) and _muon.looseId :
                has_other_muon = True
                break
        
        # Special for nanoV10 2018 Data
        if self.nanoVer == 10:
            for _ele in electrons:
                _ele_v4 = _ele.p4()
                if (_ele_v4.Pt() > 10) and (abs(_ele_v4.Eta()) < 2.5) and (_ele.mvaIso_WPL > 0.5):
                    has_ele = True
                    break
        elif (self.nanoVer == 11): # nano v11 mvaIso_Fall17V2_WPL,  in nano v13 Electron_mvaIso
            for _ele in electrons:
                _ele_v4 = _ele.p4()
                if (_ele_v4.Pt() > 10) and (abs(_ele_v4.Eta()) < 2.5) and (_ele.mvaIso > 0.5):
                    has_ele = True
                    break
        else:
            for _ele in electrons:
                _ele_v4 = _ele.p4()
                if (_ele_v4.Pt() > 10) and (abs(_ele_v4.Eta()) < 2.5) and (_ele.mvaFall17V2Iso_WPL > 0.5):
                    has_ele = True
                    break
            
        # Apply MT cut (if enabled)
        if mtCut > 0:
            if len(lt_muon_sel) >=1:
                if self.calc_MT(signalMuon_v4, self.calc_met(met)) > mtCut:
                    MT_Cut = True
        
        # Apply b tag veto (if enabled)
        if btagThreshold > 0:
            for _jet in jets:
                _j_v4 = _jet.p4()
                if (_j_v4.Pt() > 20) and (abs(_j_v4.Eta()) < 2.4) and (_jet.btagDeepFlavB > btagThreshold):
                    has_bjet = True
                    break

        # Apply MET filter
        if not flag.METFilters:
            pass_MET_filter = False

        # return result
        '''
        if (find_signal_muon) and (not has_other_muon) and (not has_ele) and (pass_MET_filter):
            if MT_Cut:
                print("This event has not pass due to MT Cut!")
                return False
            elif has_bjet:
                print("THis event has not pass due to has bjets!")
                return False
            else:
                print("This event has pass!")
                return True
        else:
            print("This event is not pass and its signal muon num is {} !".format(len(lt_muon_sel)))
            return False
        '''
        self.fill_cut('no_cut')
        if find_signal_muon:
            self.fill_cut('signal_muon')
            if not has_other_muon:
                self.fill_cut('muon_veto')
                if not has_ele:
                    self.fill_cut('ele_veto')
                    if pass_MET_filter:
                        self.fill_cut('met_filters')
                        if not MT_Cut:
                            self.fill_cut('MT_Cut')
                            if not has_bjet:
                                self.fill_cut('btag_veto')
                                return True
        

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

selection2016MC = lambda : selectionFilter(True,"2016")
selection2017MC = lambda : selectionFilter(True,"2017")
selection2018MC = lambda : selectionFilter(True,"2018")
selection2022MC = lambda : selectionFilter(True,"2022")
selection2023MC = lambda : selectionFilter(True,"2023")
selection2016data = lambda : selectionFilter(False,"2016")
selection2017data = lambda : selectionFilter(False,"2017")
selection2018data = lambda : selectionFilter(False,"2018")
selection2022data = lambda : selectionFilter(False,"2022")
selection2023data = lambda : selectionFilter(False,"2023")
