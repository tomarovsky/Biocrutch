#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
from Biocrutch.Parsers.url_parsers import sourceforge_latest_link_and_version
from Biocrutch.Routines.routine_functions import extract_file
import os

# auto for sourceforge.net
# test for samtools


def main():
    # parsing latest version and download link
    latest_tool_info = sourceforge_latest_link_and_version(args.tool)
    print('INFO:', latest_tool_info)
    tool_link = latest_tool_info[0]
    tool_archive_name = latest_tool_info[1].split()[0]
    tool_directory = tool_archive_name.rsplit('.', 2)[0]
    if '-' in tool_directory:
        version = tool_directory.split('-')[-1]
    elif '_' in tool_directory:
        version == tool_directory.split('_')[-1]
    print('LINK:', tool_link)
    print('VERSION:', tool_directory)  # to print the tool name too

    # check version by available directory
    current_directory = os.getcwd()
    os.chdir(args.working_directory + args.tool)
    print(os.getcwd())
    if os.path.isdir(tool_directory):
        print('latest tool version is installed!')
        os.chdir(current_directory)
        exit()
    else:
        print('latest tool version is NOT installed!')
        if not input('Do you want a new version? (y/n): ').lower().strip()[:1] in ['y', 'yes']:
            print('ok. mb later')
        else:
            print('module load gcc/gcc-6.3')
            os.system('module load gcc/gcc-6.3')
            print('Downloading files...')
            os.system('wget --trust-server-names ' + tool_link)
            print('The tool files have been downloaded. Unpacking ...')
            extract_file(tool_archive_name)
            print('Installation...')
            os.chdir(tool_directory)
            print(os.getcwd())
            os.system(args.install_command)
            print('Check:')
            os.system('./{} --version'.format(args.tool))
            os.chdir('..')
            os.remove(tool_archive_name)

            # create modulefile
            print('Creating a module...')
            modulefile_template = [
                '#%Module1.0#####################################################################',
                '##',
                '## modules modulefile',
                '##',
                '## modulefiles/modules.  Generated from modules.in by configure.',
                '##',
                'module-whatis   "Activates {tool} package"',
                '',
                '# for Tcl script use only',
                'module load gcc/gcc-6.3',
                'set\tversion\t{version}',
                'set\tmodroot\t/usr/share/Modules',
                'set\ttopdir\t/home/tools/{tool}/{tool_directory}',
                'prepend-path\tPATH\t$topdir',
            ]
            if not os.path.exists('/home/tools/modulefiles/{tool}'.format(tool=args.tool)):
                os.mkdir('/home/tools/modulefiles/{tool}'.format(tool=args.tool))
            with open(
                '/home/tools/modulefiles/{tool}/{tool_directory}'.format(tool=args.tool, version=version, tool_directory=tool_directory), 'w'
            ) as modulefile:
                modulefile.write('\n'.join(modulefile_template).format(tool=args.tool, version=version, tool_directory=tool_directory))
            print('Successfully installed')
            print('Use: module load {tool}/{tool_directory}'.format(tool=args.tool, tool_directory=tool_directory))

    os.chdir(current_directory)
    exit()


if __name__ == '__main__':
    parser = ArgumentParser(description='automatic parser for updating tools')

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-t', '--tool', type=str, default='samtools', help='tool available at sourceforge.net')
    group_additional.add_argument('-d', '--working-directory', type=str, default='/home/tools/', help='directory type of ./tool/tool-<version>')
    group_additional.add_argument(
        '-c',
        '--install-command',
        type=str,
        default='./configure && make',
        help="The default for installing samtools. You can use different installation scripts for specific tools. For example: 'make', 'make && make install' or others.",
    )
    args = parser.parse_args()
    main()
