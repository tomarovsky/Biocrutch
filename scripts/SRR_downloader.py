from argparse import ArgumentParser
import subprocess
from Biotoolsoup.Parsers.url_parsers import SRR_download_link


def main():
    for SRR_id in args.input:
        link = SRR_download_link(SRR_id)

        #axel - bash tool for download SRR from the link
        download_SRR_bash_command = ("axel -n 50 {}".format(link)).split()
        subprocess.call(download_SRR_bash_command)

if __name__ == "__main__":
    parser = ArgumentParser(description="")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str,
                                nargs='+', required=True,
                                help="SRR id/ids")
    args = parser.parse_args()
    main()
