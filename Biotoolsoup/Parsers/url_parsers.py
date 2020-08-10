from bs4 import BeautifulSoup
import requests
import lxml
import pandas as pd

"""Various scripts related to parsing sites"""

def SRR_download_link(SRR_id):
    """NCBI SRA downloader"""
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={}".format(SRR_id)
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    try:
        for tag_a in soup.find('tr', {'class': 'first'}).find_all('a'):
            link = tag_a['href']
            return link #first link
    except AttributeError:
        print("ERROR!!! SRR not found")

def SRR_metrics(SRR_id):
    """NCBI SRA metrics checker"""
    dictionary = {}
    url = "https://www.ncbi.nlm.nih.gov/sra/?term={}".format(SRR_id)
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    count = 0
    for SRR_id in soup.find_all('td', {'align': 'left'}):
        count += 3
        SRR_id = SRR_id.get_text()
        metrics_lst = []
        metrics = soup.find_all('td', {'align': 'right'})
        for metric in metrics:
            metrics_lst.append(metric.get_text())
        if count == 3:
            dictionary[SRR_id] = metrics_lst[:3]
        else:
            dictionary[SRR_id] = metrics_lst[count - 3 : count]
    df = pd.DataFrame(dictionary)
    return df