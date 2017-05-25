from bs4 import BeautifulSoup
import os
import re
import requests
import zipfile

datasets = [
    ('https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-84422/files/', 'E-GEOD-84422')
]

def download_file(url, folder):
    print('Downloading file: {}'.format(url))
    local_filename = folder + '/' + url.split('/')[-1]
    if not os.path.isfile(local_filename):
        r = requests.get(url, stream=True)
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
    else:
        print('  - Already exists!')
    return local_filename

def download_arrayexpress(url, dataset):
    base_url = '/'.join(url.split('/')[:3])
    html = requests.get(url)

    # Sanity check
    if html.status_code is not 200:
        return None

    soup = BeautifulSoup(html.content)

    folder = '../data/' + dataset
    if not os.path.exists(folder):
        os.makedirs(folder)

    # Get sequence data
    aefiles = [a['href'] for a in soup.find('div', {'class': 'ae-files'}).find_all('a')][1:]
    for f in aefiles:
        lf = download_file(base_url + f, folder)
        if lf[-4:] == '.zip':
            z = zipfile.ZipFile(lf, 'r')
            z.extractall(folder + '/samples')
            z.close()
            os.remove(lf)  # Remove zip file to save space

    # Get array data
    adffiles = [a['href'] for a in soup.find_all('a', href=re.compile('\\.adf.txt$'))]
    for f in adffiles:
        download_file(base_url + f, folder)

for url, dataset in datasets:
    download_arrayexpress(url, dataset)

