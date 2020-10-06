#!/usr/bin/env python3
from argparse import ArgumentParser
import subprocess
from Biocrutch.Parsers.url_parsers import SRA_download_link
from Biocrutch.Parsers.url_parsers import SRA_metrics


def main():
    for SRA_id in args.input:
        if args.download:
            link = SRA_download_link(SRA_id)
            #axel - tool for download SRA from the link
            download_SRA_bash_command = ("axel -n 50 {}".format(link)).split()
            subprocess.call(download_SRA_bash_command)

        if args.metrics:
            print("Metrics:")
            SRA_id_only = SRA_id.split("_")[0]
            metrics_frame = SRA_metrics(SRA_id_only)
            print(metrics_frame)
            print("Counting reads...")
            cat = subprocess.Popen(['cat', SRA_id], stdout=subprocess.PIPE)
            wc = subprocess.Popen(['parallel', '--pipe', 'wc', '-l'], stdin=cat.stdout, stdout=subprocess.PIPE)
            awk = subprocess.Popen(['awk', '{s+=$1} END {print s}'], stdin=wc.stdout, stdout=subprocess.PIPE)
            lines_count = awk.communicate()[0]
            # print(lines_count.decode("utf-8"))
            reads_count = int(lines_count)/4
            print("Number of reads = {}".format(reads_count))
            if int(metrics_frame[SRA_id_only][0].replace(',', '')) == int(reads_count):
                print("Yes! Reads.fastq downloaded without damage")
            else:
                print("No. Reads.fastq downloaded with damage")


if __name__ == "__main__":
    parser = ArgumentParser(description="For downloading SRA data and parsing metrics.")
    group_required = parser.add_argument_group('Options')
    group_required.add_argument('-i', '--input', type=str,
                                nargs='+', help="SRA id (or file.fastq if you want metrics)")
    group_required.add_argument('-d', '--download',
                                action="store_true",
                                help="downloader SRA id/ids")
    group_required.add_argument('-m', '--metrics',
                                action="store_true",
                                help="metrics for SRA_id.fastq and read counts for comparison")
    
    args = parser.parse_args()
    main()
