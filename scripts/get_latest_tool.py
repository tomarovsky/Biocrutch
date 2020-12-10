#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
from Biocrutch.Parsers.url_parsers import sourceforge_latest_link_and_version
import os

# auto for sourceforge.net
# test for samtools

def main():
    # parsing latest version and download link
    latest_tool_info = sourceforge_latest_link_and_version(args.input)
    print('LINK:', latest_tool_info[0])
    version = latest_tool_info[1].split()[0].rsplit('.', 2)[0] # version == dir_name
    print('VERSION:', version)
    
    # check version by availble directory
    current_directory = os.path.abspath(os.curdir)
    os.chdir(args.working_directory + args.input)
    print(os.getcwd())
    if os.path.isdir(version):
        print("latest tool version is installed!")
        os.chdir(current_directory) 
        exit()
    else:
        print("latest tool version is NOT installed!")
        answers = ['y', 'yes']
        if not input("Are you sure? (y/n): ").lower().strip()[:1] in answers:
            print('ok. mb later')
            exit()
        else:
            print("Downloading files...")
    
    os.chdir(current_directory) 

if __name__ == "__main__":
    parser = ArgumentParser(description="automatic parser for updating tools")

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-i', '--input', type=str, 
                                  default="samtools", help="tool available at sourceforge.net")
    group_additional.add_argument('-d', '--working-directory', type=str,
                                  default="/home/tools/")
    args = parser.parse_args()
    main()