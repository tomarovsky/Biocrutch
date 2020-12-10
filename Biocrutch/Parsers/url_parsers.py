__author__ = 'tomarovsky'
from bs4 import BeautifulSoup
import requests
import lxml
import pandas as pd

"""Various scripts related to parsing sites"""

def SRA_download_link(SRA_id):
    """NCBI SRA data downloader"""
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={}".format(SRA_id)
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    try:
        for tag_a in soup.find('tr', {'class': 'first'}).find_all('a'):
            link = tag_a['href']
            return link #first link
    except AttributeError:
        print("ERROR!!! SRA data not found")

def SRA_metrics(SRA_id):
    """NCBI SRA data metrics checker"""
    dictionary = {}
    url = "https://www.ncbi.nlm.nih.gov/sra/?term={}".format(SRA_id)
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    count = 0
    for SRA_id in soup.find_all('td', {'align': 'left'}):
        count += 3
        SRA_id = SRA_id.get_text()
        metrics_lst = []
        metrics = soup.find_all('td', {'align': 'right'})
        for metric in metrics:
            metrics_lst.append(metric.get_text())
        if count == 3:
            dictionary[SRA_id] = metrics_lst[:3]
        else:
            dictionary[SRA_id] = metrics_lst[count - 3 : count]
    df = pd.DataFrame(dictionary)
    return df

def sourceforge_latest_link_and_version(tool_name):
    """get version number and download link"""
    url = "https://sourceforge.net/projects/{}/files/".format(tool_name)
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    files = soup.find('div', {'id': "files"})
    sam = files.find('div', {'class' : "btn-set"}).a
    link = "https://sourceforge.net" + sam.get('href')
    # print('download link:', link)
    print('ACTUAL INFO: ', sam.get('title').split('/')[-1])
    version_info = sam.find('span', {'class' : 'sub-label'}).get_text()
    # print('latest tool version:', version_info)
    return (link, version_info)