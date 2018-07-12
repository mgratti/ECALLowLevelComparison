
from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph, TH1
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta, kTRUE, kFALSE
from ROOT import TEfficiency
import itertools
import sys
import os
sys.path.append('{user}/plotting/myplotting'.format(user=os.environ['HOME']))
from spares import *
import graphUtils as gU


class EoverEtrueAnalysisResult(object):
  def __init__(self, det=0, eff=0, sigma=0, rms=0, iseeding=0, igathering=0):
    self.det = det
    self.eff = eff
    self.sigma = sigma
    self.rms = rms
    self.iseeding = iseeding
    self.igathering = igathering


def makeEoverEtrueAnalysis(inputfile, det, iseeding, igathering, nevts, outputdir, outfile, fitoutfile):

  result = EoverEtrueAnalysisResult()
  if not os.path.isfile(inputfile):
    return False,result



  ###################
  # READ
  ##################
  histoname = 'h_PFclusters_genMatched_{det}_eOverEtrue'.format(det=det)
  inputdir = 'ecalnoisestudy'
  subdir = 'PFClusters'
  print inputfile
  f=TFile(inputfile, 'READ')
  histo=f.Get('{}/{}/{}'.format(inputdir,subdir,histoname))
  if not histo: return False,result

  histo.SetMarkerStyle(20)
  histo.GetXaxis().SetTitle('E_{{PFcluster}} / E_{{True}}, {det} (GeV)'.format(det=det))
  histo.GetYaxis().SetTitle('Entries')
  histo.GetXaxis().SetRangeUser(0, 2)
  #histo.GetYaxis().SetRangeUser(0., 1600)
  #if 'EE' in det:
  #  histo.GetYaxis().SetRangeUser(0., 600)

  # better to avoid setting the range, since the mean calculation changes
  #histo.GetXaxis().SetRangeUser(xrange[0], xrange[1])

  ###################
  # FIT
  ##################
  f1 = TF1('f1','crystalball',0.4, 2.)
  f1.SetParameters(200, 1, 0.05, 3, 2) # my guess: constant (normalization)=integral, mean = 1, sigma = 0.1, alpha (quanto lontano dal picco si innesta la coda) = 0.7, N = 0.5 (lunghezza della coda(?)
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
  #fit_params = [ ('Param {}'.format(i),'{:.2f}'.format(f1.GetParameter(i)), '{:.2f}'.format(f1.GetParError(i)) ) for i in range(0,5)]
  #ffitout = open(det + '_' + fitoutfile , 'a')
  #ffitout.write('\n\nFit results for seeding={} gathering={} subdet={}:\n'.format(iseeding, igathering, det))
  #ffitout.write('\nChi2/Ndf=' + str(f1.GetChisquare()) + '/' + str(f1.GetNDF()) + '\n')
  #par_string = '\n'.join("%s: val=%s  err=%s" % tup for tup in fit_params)
  #ffitout.write(par_string)

  ###################
  # Efficiency
  ##################
  hpass = TH1F('hpass', 'hpass', 1, 0., 1.)
  Npass = histo.GetEntries()
  for i in range(0,int(Npass)):
    hpass.Fill(0.5)

  htot = f.Get('{}/{}/{}'.format(inputdir,'general', 'h_genP_n{d}'.format(d=det)))

  print 'eff calculated from npass {} over ntot {}'.format(hpass.GetEntries(),htot.GetEntries())

  if TEfficiency.CheckConsistency(hpass, htot):
    pEff = TEfficiency(hpass, htot) # default stat option is clopper pearson
    eff=pEff.GetEfficiency(1)
    erru=pEff.GetEfficiencyErrorUp(1)
    errd=pEff.GetEfficiencyErrorLow(1)
  else:
    eff = 1.0
    erru = 1.0
    errd = 1.0
  eff_label = '{:.3f}+{:.3f}-{:.3f}'.format(eff,erru,errd)
  #fout = open(det + '_' + outfile, 'a')
  #fout.write('Total efficiency for seeding={} gathering={}:  {} \n'.format(iseeding,igathering,eff_label))
  #fout.close()

  eff_label = 'Eff={}'.format(eff_label)
  defaultLabels([eff_label], x=0.62, y=0.65, spacing = 0.04, size = 0.06, dx = 0.12)


  ###################
  # RESULTS
  ##################
  c.SaveAs('{o}/EoverEtrue_{det}_{s}.pdf'.format(o=outputdir,det=det, s=suffix))
  c.SaveAs('{o}/EoverEtrue_{det}_{s}.png'.format(o=outputdir,det=det, s=suffix))

  return True,EoverEtrueAnalysisResult(det=det, eff=eff, sigma=f1.GetParameter(2), rms=histo.GetRMS(), iseeding=iseeding, igathering=igathering)


if __name__ == "__main__":

  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  #gStyle.SetOptStat('emMrRo')
  TH1.StatOverflows(kTRUE) # if false stats will be calculated without overflow, must be set also at filling time


  ####################################
  ## Define input and output
  ####################################

  version = 'vprodV1_ecalV5'
  inputfile = '../test/outputfiles/test_photonGun_seed{s}_gather{g}_{v}_numEvent{n}.root'

  params = {}
  params["nevts"] =     [10000]
  params["gathering"] = [1.0, 2.0, 5.0, 10., 0.5, 0.2, 0.1]# [1.0, 2.0, 5.0, 10., 0.5, 0.2, 0.1] # multiplier
  params["seeding"] =   [1.0, 2.0, 5.0, 10., 0.5, 0.2, 0.1] # multiplier#
  parameters_set = list(itertools.product(params["nevts"],params["seeding"],params["gathering"] ))

  results = []
  outputdir = 'plots/anaEoverEtrue_{v}'.format(v=version)

  os.system('mkdir {}'.format(outputdir))

  ####################################
  ## Do the analysis for each seeding-gathering threshold
  ####################################
  for iset in parameters_set:
    inevts,iseeding,igathering = iset
    inputfilename = inputfile.format(s=iseeding, g=igathering, n=inevts, v=version)
    for det in ['EB']:
      ret,result = makeEoverEtrueAnalysis(inputfilename, det, iseeding, igathering, inevts, outputdir, outfile='Efficiency_EoverEtrue.txt', fitoutfile='fit_EoverEtrue.txt')
      if ret:
        results.append(result)
      else:
        print 'No result for' , iset, '  skipping'

  ###################################
  ## Plot the results
  ##################################
  effs = []
  resofit = []
  resorms = []
  for result in results:
    if result.det != 'EB': continue
    effs.append(result.eff)
    resofit.append(result.sigma)
    resorms.append(result.rms)

  if len(effs)!=0:

    g1 = gU.makeGraph( xs=effs, ys=resofit, xtitle='N_{reco}/N_{gen}', ytitle='#sigma(E) (Fit)', label='Et true = 9-10 GeV', color=kPink, style=20)
    gU.makePlot(graphs=[g1], plotName='EffVsReso', outputdir=outputdir)

    g2 = gU.makeGraph( xs=effs, ys=resorms, xtitle='N_{reco}/N_{gen}', ytitle='#sigma(E) (RMS)', label='Et true = 9-10 GeV', color=kPink, style=20)
    gU.makePlot(graphs=[g2], plotName='EffVsResoRMS', outputdir=outputdir)
