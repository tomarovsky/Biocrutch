from bs4 import BeautifulSoup
import requests
import lxml

"""Various scripts related to parsing sites"""

def SRR_download_link(SRR_id):
    """NCBI SRA"""
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={}".format(SRR_id)
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    try:
        for tag_a in soup.find('tr', {'class': 'first'}).find_all('a'):
            link = tag_a['href']
            return link #first link
    except AttributeError:
        print("ERROR!!! SRR not found")