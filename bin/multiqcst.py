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
def ck_filegeneration(inputpath,filenumber):
    if len(os.listdir(inputpath)) >= int(filenumber):
        #print("complete")
        complete =True
        return complete
    else:
        #print("Still waiting for fastqc")
        complete =False
        return complete


#The waiting function
def waiting_fastqc(inputpath,multiqcpath,filenumber):
    complete=False
    while not complete :
        complete=ck_filegeneration(inputpath,filenumber)
        if complete == False:
            print("Waiting for fastqc to finish .....")
            time.sleep(10)
            return complete
        else:
            complete=ck_filegeneration(inputpath,filenumber)
            return complete


#The multiqc after fastqc function
def launch_multiqc(inputpath,multiqcpath):
    os.system("cd " +  multiqcpath)
    os.system("multiqc . ")


#The check
def check_fastqc_finish():
    #Arg call
    parser = argparse.ArgumentParser()
    parser.add_argument("inputpath", type=str, help="Set fastqc Path.")
    parser.add_argument("multiqcpath", type=str, help="Set MultiQC Path.")
    parser.add_argument("filenumber", type=int, help="Number of file to wait for.")
    parser.parse_args()
    if (len(sys.argv) < 4):
        print("You have ommitted an argument")
        print ("Please put the right arguments and try again :)")
        sys.exit('Usage: python %s --help' % sys.argv[0])

    elif (len(sys.argv) == 4):
        inputpath = sys.argv[1]
        multiqcpath = sys.argv[2]
        filenumber = sys.argv[3]

        # Call the waiting function
        complete= waiting_fastqc(inputpath,multiqcpath,filenumber)
        if complete == True:
            #call the multiqc after trimming function
            launch_multiqc(inputpath,multiqcpath)

    else:
        #print (indir)
        print("Something went wrong please contact the developer ...")
        sys.exit()

#Call function
check_fastqc_finish()