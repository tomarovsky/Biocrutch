import os
import argparse
from Pipelines import FilteringPipeline
from RouToolPa.Collections.General import IdList
from RouToolPa.Routines.File import check_path


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with samples")
parser.add_argument("-x", "--general_stat_file", action="store", dest="general_stat_file", required=True,
                    help="File to write general statistics about filtration")
parser.add_argument("-s", "--samples", action="store", dest="samples", type=lambda s: s.split(","),
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples")
parser.add_argument("-f", "--sample_file", action="store", dest="sample_file",
                    help="File with sample ids to handle. No effect if -s/--samples option is set")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use. Default - 1.")

parser.add_argument("-a", "--adapters", action="store", dest="adapters", type=os.path.abspath,
                    help="File with adapters to trim by Trimmomatic")
parser.add_argument("-k", "--adapter_kmers", action="store", dest="adapter_kmers", type=os.path.abspath,
                    help="File with adapter k-mers for Coockiecutter")

parser.add_argument("-m", "--mismatch_number", action="store", dest="mismatch_number", type=int, default=2,
                    help="Number of mismatches in adapter seed. Works only if -a/--adapters option is set. Default - 2.")
parser.add_argument("-p", "--pe_score", action="store", dest="pe_score", type=int, default=30,
                    help="PE reads adapter score. Works only if -a/--adapters option is set. Default - 30.")
parser.add_argument("-e", "--se_score", action="store", dest="se_score", type=int, default=10,
                    help="SE reads adapter score. Works only if -a/--adapters option is set. Default - 10.")
parser.add_argument("-n", "--min_adapter_len", action="store", dest="min_adapter_len", type=int, default=1,
                    help="Minimum length of adapter fragment. Works only if -a/--adapters option is set. Default - 1.")
parser.add_argument("-g", "--sliding_window_size", action="store", dest="sliding_window_size", type=int,
                    help="Size of sliding window when checking quality. "
                         "If not set - reads will be filtered by mean quality")
parser.add_argument("-q", "--average_quality_threshold", action="store", dest="average_quality_threshold", default=20,
                    type=int,
                    help="Quality threshold for sliding window or whole read."
                         "Depends on -q/--average_quality_threshold option.Default - 20.")
parser.add_argument("-b", "--base_quality", action="store", dest="base_quality", default="phred33",
                    help="Type of base quality. Possible variants: phred33, phred64. Default - phred33 ")
parser.add_argument("-l", "--min_length", action="store", dest="min_len", type=int, default=50,
                    help="Minimum length of read to retain. Default - 50")
parser.add_argument("-j", "--trimmomatic_dir", action="store", dest="trimmomatic_dir", default="",
                    help="Path to Trimmomatic directory")
parser.add_argument("-c", "--trimmer_dir", action="store", dest="trimmer_dir", default="",
                    help="Path to Trimmer directory")

parser.add_argument("-r", "--remove_intermediate_files", action="store_true",
                    dest="remove_intermediate_files", default=False,
                    help="Remove intermediate files")
parser.add_argument("--stat_off", action="store_true",
                    dest="stat_off", default=False,
                    help="Turn off calculation of statistics. Default: False")

args = parser.parse_args()

FilteringPipeline.stirka_trimmomatic(args.samples_dir, args.output_dir, args.adapter_kmers, args.adapters,
                                     args.general_stat_file,
                                     samples_to_handle=args.samples if args.samples else IdList(filename=args.sample_file) if args.sample_file else None,
                                     threads=args.threads, trimmomatic_dir=args.trimmomatic_dir,
                                     trimmer_dir=args.trimmer_dir,
                                     mismatch_number=args.mismatch_number, pe_reads_score=args.pe_score,
                                     se_read_score=args.se_score, min_adapter_len=args.min_adapter_len,
                                     sliding_window_size=args.sliding_window_size,
                                     average_quality_threshold=args.average_quality_threshold,
                                     leading_base_quality_threshold=None, trailing_base_quality_threshold=None,
                                     crop_length=None, head_crop_length=None, min_len=args.min_len,
                                     base_quality=args.base_quality,
                                     remove_intermediate_files=args.remove_intermediate_files,
                                     stat_off=args.stat_off)