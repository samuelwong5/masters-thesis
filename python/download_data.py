from bs4 import BeautifulSoup
import os
import re
import requests
import sys
import zipfile


datasets = [
    ('https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-84422/files/', 'E-GEOD-84422'),
    ('https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-63063/files/', 'E-GEOD-63063'),
    ('https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-48350/files/', 'E-GEOD-48350'),
    ('https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-28894/files/', 'E-GEOD-28894'),
    ('https://www.ebi.ac.uk/arrayexpress/experiments/E-AFMX-6/files/', 'E-AFMX-6'),
    ('https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-3790/files/', 'E-GEOD-3790')
]

def print_nonewline(s):
    sys.stdout.write(s)
    sys.stdout.flush()

def download_file(url, folder):
    print('Downloading file: {}'.format(url.split('/')[-1]))
    local_filename = folder + '/' + url.split('/')[-1]
    if not os.path.isfile(local_filename):
        r = requests.get(url, stream=True)
        size = int(r.headers['Content-length']) // 1024 // 1024
        with open(local_filename, 'wb') as f:
            counter = 0
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    counter += 1
                    print_nonewline('\r  - {}/{}MB'.format(counter // 1024, size))
            print()
    else:
        print('  - File already exists, skipping.')
        return None
    return local_filename


def download_arrayexpress(url, dataset):
    base_url = '/'.join(url.split('/')[:3])
    html = requests.get(url)

    # Sanity check
    if html.status_code is not 200:
        return None

    soup = BeautifulSoup(html.content, 'html5lib')

    folder = '../data/' + dataset
    if not os.path.exists(folder):
        os.makedirs(folder)
    print('Downloading dataset {} to {}'.format(dataset, folder))

    # Get sequence data
    aefiles = [a['href'] for a in soup.find('div', {'class': 'ae-files'}).find_all('a')]
    aefiles = list(filter(lambda x: not 'raw' in x, aefiles))
   
    for f in aefiles:
        lf = download_file(base_url + f, folder)
        if lf and lf[-4:] == '.zip':
            z = zipfile.ZipFile(lf, 'r')
            z.extractall(folder + '/samples')
            z.close()
            # os.remove(lf)  # Remove zip file to save space

    # Get array data
    adffiles = [a['href'] for a in soup.find_all('a', href=re.compile('\\.adf.txt$'))]
    for f in adffiles:
        lf = download_file(base_url + f, folder)
        if lf:
            truncate_adf(lf)


def truncate_adf(adf_file):
    try:
        f = open(adf_file, 'r+')
        d = f.readlines()
        f.seek(0)
        index = d.index('[main]\n')
        print(adf_file, index)
        for i in d[index+1:]:
            f.write(i)
            f.truncate()
    except ValueError:
        pass # ADF file does not have metadata
    finally:
        f.close()


if __name__ == '__main__':
    for url, dataset in datasets:
        download_arrayexpress(url, dataset)
