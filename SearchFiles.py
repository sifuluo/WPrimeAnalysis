import os

samplebasepaths = '/fdata/hepx/store/user/aoverton0342/madGraph/ak4/'
folders = ['TDual_FormerLeptonic','TDual_LatterLeptonic']
for folder in folders:
  f = open(folder+".txt","w+")
  basepath = samplebasepaths + folder + "/Events/"
  irun = 1
  while os.path.exists(basepath + "run_{:0>2d}/".format(irun)):
    runpath = basepath + "run_{:0>2d}/".format(irun)
    for itag in range(1,4):
      file = runpath + "tag_{}_delphes_events.root".format(itag)
      if os.path.exists(file):
        f.write(file+"\n")
    irun+=1
