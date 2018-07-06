
from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph, TH1
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta, kTRUE, kFALSE
from ROOT import TEfficiency
import itertools
import sys
import os
sys.path.append('{user}/plotting/myplotting'.format(user=os.environ['HOME']))
from spares import *

def makeEoverEtrueAnalysis(inputfile, iseeding, igathering, nevts, outfile, fitoutfile):


  inputfile = inputfile.format(s=iseeding, g=igathering, n=inevts, v=version)
  for det in ['EB', 'EEP', 'EEM']:

    histoname = 'h_PFclusters_genMatched_{det}_eOverEtrue'.format(det=det)
    inputdir = 'ecalnoisestudy'
    subdir = 'PFClusters'
    print inputfile
    f=TFile(inputfile, 'READ')
    histo=f.Get('{}/{}/{}'.format(inputdir,subdir,histoname))

    histo.SetMarkerStyle(20)
    histo.GetXaxis().SetTitle('E_{{PFcluster}} / E_{{True}}, {det} (GeV)'.format(det=det))
    histo.GetYaxis().SetTitle('Entries')
    histo.GetXaxis().SetRangeUser(0, 2)
    histo.GetYaxis().SetRangeUser(0., 1600)
    if 'EE' in det:
      histo.GetYaxis().SetRangeUser(0., 600)

    # better to avoid setting the range, since the mean calculation changes
    #histo.GetXaxis().SetRangeUser(xrange[0], xrange[1])

    # fit
    f1 = TF1('f1','crystalball',0.4, 2.)
    f1.SetParameters(60, 1, 0.1, 0.7, 10) # my guess: constant= 60, mean = 1, sigma = 0.1, alpha = 0.7, N = 0.5
    f1.SetLineColor(kRed)

    c = TCanvas()
    histo.Draw('PE')
    #fitresult = histo.Fit(f1, 'RS')
    fitresult = histo.Fit(f1, 'R') # L for loglikelihood ,
    f1.Draw('same')
    suffix = 'seed{s}_gather{g}'.format(s=iseeding,g=igathering)
    # save later

    #fo.cd()
    #fitresult.Write()
    #fo.Close()

    # Get the fitted function parameters and write them to txt file
    fit_params = [ ('Param {}'.format(i),'{:.2f}'.format(f1.GetParameter(i)), '{:.2f}'.format(f1.GetParError(i)) ) for i in range(0,5)]
    ffitout = open(det + '_' + fitoutfile , 'a')
    ffitout.write('\n\nFit results for seeding={} gathering={} subdet={}:\n'.format(iseeding, igathering, det))
    ffitout.write('\nChi2/Ndf=' + str(f1.GetChisquare()) + '/' + str(f1.GetNDF()) + '\n')
    par_string = '\n'.join("%s: val=%s  err=%s" % tup for tup in fit_params)
    ffitout.write(par_string)

    # efficiency stuff
    hpass = TH1F('hpass', 'hpass', 1, 0., 1.)
    Npass = histo.GetEntries()
    for i in range(0,int(Npass)):
      hpass.Fill(0.5)

    htot = f.Get('{}/{}/{}'.format(inputdir,'general', 'h_genP_n{d}'.format(d=det)))
    #print 'bin contents', hpass.GetBinContent(1), htot.GetBinContent(1),
    fout = open(det + '_' + outfile, 'a')
    if TEfficiency.CheckConsistency(hpass, htot):
      pEff = TEfficiency(hpass, htot) # default stat option is clopper pearson
      eff=pEff.GetEfficiency(1)
      erru=pEff.GetEfficiencyErrorUp(1)
      errd=pEff.GetEfficiencyErrorLow(1)
      eff_label = '{:.3f}+{:.3f}-{:.3f}'.format(eff,erru,errd)
      fout.write('Total efficiency for seeding={} gathering={}:  {} \n'.format(iseeding,igathering,eff_label))
    fout.close()

    eff_label = 'Eff={}'.format(eff_label)
    defaultLabels([eff_label], x=0.62, y=0.65, spacing = 0.04, size = 0.06, dx = 0.12)

    c.SaveAs('plots/EoverEtrue_{det}_{s}.pdf'.format(det=det, s=suffix))
    c.SaveAs('plots/EoverEtrue_{det}_{s}.png'.format(det=det, s=suffix))


if __name__ == "__main__":

  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  #gStyle.SetOptStat('emMrRo')
  TH1.StatOverflows(kTRUE) # if false stats will be calculated without overflow, must be set also at filling time

  outfile = 'Efficiency_EoverEtrue.txt'
  fitoutfile = 'fit_EoverEtrue.txt'
  version = 'v2'

  #fo=TFile(rootoutfile, 'NEW')

  params = {}
  params["nevts"] =     [10000]
  params["gathering"] = [0.5,1.0, 2.0]#[1.0] # multiplier
  params["seeding"] =   [0.5,1.0, 2.0] # multiplier
  parameters_set = list(itertools.product(params["nevts"],params["seeding"],params["gathering"] ))

  for iset in parameters_set:
    inevts,iseeding,igathering= iset
    inputfile = '../test/outputfiles/test_photonGun_seed{s}_GATHER{g}_{v}_numEvent{n}.root'
    makeEoverEtrueAnalysis(inputfile, iseeding, igathering, inevts, outfile, fitoutfile)
