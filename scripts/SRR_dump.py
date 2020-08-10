from argparse import ArgumentParser
import subprocess
from Biotoolsoup.Parsers.url_parsers import SRR_download_link
from Biotoolsoup.Parsers.url_parsers import SRR_metrics


def main():
    for SRR_id in args.input:
        if args.download:
            link = SRR_download_link(SRR_id)
            #axel - tool for download SRR from the link
            download_SRR_bash_command = ("axel -n 50 {}".format(link)).split()
            subprocess.call(download_SRR_bash_command)

        if args.metrics:
            print("Metrics:")
            metrics_frame = SRR_metrics(SRR_id)
            print(metrics_frame)
            print("Counting reads...")
            reads_count = ("wc -l {} | awk '{$1=$1/4; print}'".format(SRR_id)).split()
            print("Number of reads = {}".format(reads_count))
            if int(metrics_frame[SRR_id][0].replace(',', '')) == int(reads_count):
                print("Yes! Reads.fastq downloaded without damage")


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
