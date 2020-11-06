#!/usr/bin/python
"""
Script for path checking
By Najlabioinfo Oct20
"""

import argparse
import sys
import os
#import pathlib

def isitdir(directory):
    if (os.path.isdir(directory)) == True:
        return  True
    else:
        print("Please check your path")
        return False

def createdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print("Your path dosen t exist , it will be created")
    else:
        print("OK: Your is correct.")

def check_it():
    #""" Check the arg from nextflow script """


    #Arg call
    parser = argparse.ArgumentParser()
    parser.add_argument("username", type=str, help="Script user.")
    parser.add_argument("indir", type=str, help="input dir with file to check.")
    parser.add_argument("outdir", type=str, help="output dir for report edition.")
    parser.parse_args()
    if (len(sys.argv) < 4):
        print("you have ommitted an argument")
        print ("Please put the right arguments and try again :)")
        sys.exit('Usage: python %s --help' % sys.argv[0])

    elif (len(sys.argv) == 4):
        username = sys.argv[1]
        indir = sys.argv[2]
        outdir = sys.argv[3]

        print("\n************************************************")
        print("Hello " + username + " from NajlaBioinfo !! ")
        print("************************************************")
        
        #Input dir
        print("InputDir : " + indir)
        if isitdir(indir)==True:
            print("OK: Correct Path")
            if len(os.listdir(indir))== 0:
                print("Your directory is empty put some data and try again :) !!")
            else:
                
                print("OK: " + str(len([f for f in os.listdir(indir)if os.path.isfile(os.path.join(indir, f))])) + " files ready for chek ;) !!")

        else:
            #print (indir)
            print("Please check your path")
            sys.exit()
        #Output dir
        print("OutputDir: " + outdir)
        if isitdir(outdir) == False:
            #print(outdir)
            createdir(outdir)
        else:
            print("OK: Correct Path")
        print("************************************************\n")
    
    else:
        print ("Please put the right arguments and try again :)")
        sys.exit('Usage: python %s --help' % sys.argv[0])
    
    #Footer
    #print("--------------------------------------------------------------------------------------------")
    #print("Thank you for using  argument check script \nMore infos : bhndevtools@gmail.com")
    #print("--------------------------------------------------------------------------------------------")


#Launch check arguments
check_it()