#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
from Biocrutch.Parsers.url_parsers import sourceforge_latest_link_and_version
import os

# auto for sourceforge.net
# test for samtools

def main():
    current_directory = os.path.abspath(os.curdir)
    os.chdir(args.working_directory + args.input)
    print(os.getcwd())
    
    #parsing latest version and download link
    latest_tool_info = sourceforge_latest_link_and_version(args.input)
    print('download link:', latest_tool_info[0])
    print('latest tool version info:', latest_tool_info[1])
    version = latest_tool_info[1].split()[0].rsplit('.', 2)[0]
    print('latest tool version', version)
    
    # check version
    #нужно делать на серве
    
    
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