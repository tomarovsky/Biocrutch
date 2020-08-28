from argparse import ArgumentParser
import subprocess
from Biocrutch.Parsers.url_parsers import SRR_download_link
from Biocrutch.Parsers.url_parsers import SRR_metrics


def main():
    for SRR_id in args.input:
        if args.download:
            link = SRR_download_link(SRR_id)
            #axel - tool for download SRR from the link
            download_SRR_bash_command = ("axel -n 50 {}".format(link)).split()
            subprocess.call(download_SRR_bash_command)

        if args.metrics:
            print("Metrics:")
            SRR_id_only = SRR_id.split("_")[0]
            metrics_frame = SRR_metrics(SRR_id_only)
            print(metrics_frame)
            print("Counting reads...")
            cat = subprocess.Popen(['cat', SRR_id], stdout=subprocess.PIPE)
            wc = subprocess.Popen(['parallel', '--pipe', 'wc', '-l'], stdin=cat.stdout, stdout=subprocess.PIPE)
            awk = subprocess.Popen(['awk', '{s+=$1} END {print s}'], stdin=wc.stdout, stdout=subprocess.PIPE)
            lines_count = awk.communicate()[0]
            # print(lines_count.decode("utf-8"))
            reads_count = int(lines_count)/4
            print("Number of reads = {}".format(reads_count))
            if int(metrics_frame[SRR_id_only][0].replace(',', '')) == int(reads_count):
                print("Yes! Reads.fastq downloaded without damage")
            else:
                print("No. Reads.fastq downloaded with damage")


if __name__ == "__main__":
    parser = ArgumentParser(description="For downloading files.sra, parsing metrics.")
    group_required = parser.add_argument_group('Options')
    group_required.add_argument('-i', '--input', type=str,
                                nargs='+', help="SRR id (or file.fastq if you want metrics)")
    group_required.add_argument('-d', '--download',
                                action="store_true",
                                help="downloader SRR id/ids")
    group_required.add_argument('-m', '--metrics',
                                action="store_true",
                                help="metrics for SRR_id.fastq and read counts for comparison")
    
    args = parser.parse_args()
    main()
