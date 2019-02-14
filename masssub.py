import os

basepath = "/fdata/hepx/store/user/aoverton0342/madGraph/ak4/"
# folders = "TDual_FormerLeptonic/"
Sampletype = 2
# 0 for TDual_FormerLeptonic
# 1 for TDual_LatterLeptonic
# 2 for Background
Samplename = ""
if Sampletype == 0:
  Samplename = "TDual_FormerLeptonic"
elif Sampletype == 1:
  Samplename = "TDual_LatterLeptonic"
else:
  Samplename = "ttbar"

basepath = basepath + Samplename + "/Events/"

templatef = open("submits/temp.slurm","r")
templatec = templatef.readlines()
templatef.close()

irun = 0
while True:
  irun +=1
  folder_path = basepath + ("run_%.2d/" % irun)
  if not os.path.isdir(folder_path):
    break
maxrun = irun - 1
command1 = "#SBATCH -a 0-%d"% (maxrun)
command2 = "root -l -b \"wprime.cc+(%d,$SLURM_ARRAY_TASK_ID)\" \n" % (Sampletype)
command3 = "#SBATCH -o slurmlog/%s_log_%%A-%%a.out \n" % (Samplename)
command4 = "#SBATCH -J WPrime_%s" % (Samplename)
contents = templatec
contents.insert(16, command2)
contents.insert(7, command3)
contents.insert(6, command1)
contents.insert(1,command4)
contents = "".join(contents)
batchfile = open("submit_%i.slurm" % (Sampletype),"w+")
batchfile.write(contents)
batchfile.close()
  # for itag in range(1,3):
  #   file_path = folder_path + ("tag_%.1i_delphes_events.root" % itag)
  #   if (os.path.isfile(file_path)):
  #     command = "root -l -b \"wprime.cc+(%d,%d,%d)\"" % (Sampletype, irun, itag)
  #     contents = templatec
  #     contents.insert(16, command)
  #     contents = "".join(contents)
  #
  #     print(file_path)
  #     batchfile = open("submits/submit_%i.slurm" % (irun),"w+")
  #     batchfile.write(contents)
  #     batchfile.close()
