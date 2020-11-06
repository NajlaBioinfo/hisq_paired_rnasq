#!/usr/bin/python
"""
Script for fastqc after trimming
Check if cutadapt tasks were complete 
or Waiting for them
By Najlabioinfo Oct20
"""

import argparse
import sys
import os
import time

#The file generatiion function
def ck_filegeneration(trimmingpath):
    if len(os.listdir(trimmingpath)) >= 6:
        #print("complete")
        complete =True
        return complete
    else:
        #print("Still waiting for cutadapt")
        complete =False
        return complete


#The waiting function
def waiting_cutadapt(trimmingpath,multiqcpath):
    complete=False
    while not complete :
        complete=ck_filegeneration(trimmingpath)
        if complete == False:
            print("Waiting for cutadapt to finish .....")
            time.sleep(10)
            return complete
        else:
            complete=ck_filegeneration(trimmingpath)
            return complete


#The fastqc after trimming function
def launch_fastqat(trimmingpath,multiqcpath):
    os.system("cd " +  trimmingpath)
    os.system("fastqc  * ")
    os.system("mv *fastqc* " +  multiqcpath)


#The check
def check_cutadapt_finish():
    #Arg call
    parser = argparse.ArgumentParser()
    parser.add_argument("trimmingpath", type=str, help="Set Trimming Path.")
    parser.add_argument("multiqcpath", type=str, help="Set MultiQC Path")
    parser.parse_args()
    if (len(sys.argv) < 3):
        print("you have ommitted an argument")
        print ("Please put the right arguments and try again :)")
        sys.exit('Usage: python %s --help' % sys.argv[0])

    elif (len(sys.argv) == 3):
        trimmingpath = sys.argv[1]
        multiqcpath = sys.argv[2]

        # Call the waiting function
        complete= waiting_cutadapt(trimmingpath,multiqcpath)
        if complete == True:
            #call the fastqc after trimming function
            launch_fastqat(trimmingpath,multiqcpath)

    else:
        #print (indir)
        print("Something went wrong please contact the developer ...")
        sys.exit()

#Call function
check_cutadapt_finish()