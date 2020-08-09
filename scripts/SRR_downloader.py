from argparse import ArgumentParser
from bs4 import BeautifulSoup
import requests
import lxml
import subprocess


def sra_link_parser(url):
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    for tag_a in soup.find('tr', {'class': 'first'}).find_all('a'):
        link = tag_a['href']
        return link


def main():
    for sra_id in args.input:
        url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={}".format(sra_id)
        link = sra_link_parser(url)
        download_sra_bash_command = "axel -n 50 {}".format(link)
        subprocess.run(download_sra_bash_command)


if __name__ == "__main__":
    parser = ArgumentParser(description="")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str,
                                nargs='+', required=True, help="SRR id")
    args = parser.parse_args()
    main()
