from argparse import ArgumentParser
import subprocess
from Biotoolsoup.Parsers.url_parsers import SRR_download_link
from Biotoolsoup.Parsers.url_parsers import SRR_metrics
import os

def main():
    for SRR_id in args.input:
        if args.download:
            link = SRR_download_link(SRR_id)
            #axel - tool for download SRR from the link
            download_SRR_bash_command = ("axel -n 50 {}".format(link)).split()
            subprocess.call(download_SRR_bash_command)

        if args.metrics:
            print("Metrics:")
            metrics_frame = SRR_metrics(SRR_id.split("_")[0])
            print(metrics_frame)
            print("Counting reads...")
            lines_count_command = ("wc -l {}".format(SRR_id)).split()
            lines_count_run = subprocess.run(lines_count_command,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE,
                                             universal_newlines=True)
            lines_count = lines_count_run.stdout.split()[0]
            print("Number of lines = {}".format(lines_count))
            if (4 * int(metrics_frame[SRR_id][0].replace(',', ''))) == int(lines_count):
                print("Yes! Reads.fastq downloaded without damage")
            else:
                print("No. Reads.fastq downloaded with damage")


if __name__ == "__main__":
    parser = ArgumentParser(description="for parsing links and download links using the axel.")
    group_required = parser.add_argument_group('Options')
    group_required.add_argument('-i', '--input', type=str,
                                nargs='+', help="SRR id or file.fastq if you want metrics")
    group_required.add_argument('-d', '--download',
                                action="store_true",
                                help="downloader SRR id/ids")
    group_required.add_argument('-m', '--metrics',
                                action="store_true",
                                help="metrics for SRR_id.fastq")
    
    args = parser.parse_args()
    main()
