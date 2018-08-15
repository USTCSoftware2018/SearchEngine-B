#!/usr/bin/env python3

from typing import List, Dict
import requests
import json
import csv
from bs4 import BeautifulSoup
from urllib.parse import urlencode, parse_qs


class Result:
    def __init__(self, title, url, summary):
        self.title = title
        self.url = url
        self.summary = summary


class Source:
    """
    An abstract base class for all data sources.

    Required methods are `count` and `search`.
    """
    def count(self, keyword: str) -> int:
        """
        Count the number of all possible results.

        :return: Total number of results.
        """
        raise NotImplementedError()

    def search(self, keyword: str, rng: range) -> List[Result]:
        """
        Search a certain data source for a specific keyword within some range.

        The range should fall in [0, count()).

        :return: List of `Result`s.
        """
        raise NotImplementedError()


class HumanGeneSource(Source):
    """
    Source: https://ghr.nlm.nih.gov/gene
    """
    def __init__(self):
        self._cached_first_soup = None   # Cached BeautifulSoup for the first page of results

    def soupize(self, keyword: str, params: Dict[str, str] = {}) -> BeautifulSoup:
        query_string = urlencode({'query': keyword, **params})
        query_url = 'https://ghr.nlm.nih.gov/search?%s' % query_string
        response = requests.get(query_url)
        return BeautifulSoup(response.content, features='html.parser')

    def search_oneshot(self, keyword: str, start: int) -> List[Result]:
        """
        One-time search, starting from @start, for 10 items.
        """
        if start == 0 and self._cached_first_soup:
            soup = self._cached_first_soup
        else:
            soup = self.soupize(keyword, {'start': start})

        ul = soup.find(class_='search-results')
        if not ul:
            return []

        results = []
        for li in ul.children:
            result = Result(
                li.a.text,
                'https://ghr.nlm.nih.gov{}'.format(li.a['href']),
                li.find(class_='sample-content').text
            )
            results.append(result)

        return results

    def count(self, keyword: str) -> int:
        soup = self.soupize(keyword)
        self._cached_first_soup = soup
        return int(soup.find(class_='ss-item').text[5:-1])

    def search(self, keyword: str, rng: range) -> List[Result]:
        results = []
        offset = rng.start
        step = 10

        while offset < rng.stop:
            oneshot = self.search_oneshot(keyword, offset)
            results.extend(oneshot[:min(rng.stop-offset, step)])
            if len(oneshot) < step:
                break
            offset += len(oneshot)

        return results


class iGEMPartSource(Source):
    """
    Source: http://parts.igem.org/Main_Page
    """
    pass


class RCSBSource(Source):
    """
    Source: http://www.rcsb.org/
    """
    def __init__(self):
        self.cached_ids = {}

    def query_all(self, keyword: str) -> List[str]:
        if keyword in self.cached_ids:
            return self.cached_ids[keyword]

        data = '<orgPdbQuery><queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType><keywords>%s</keywords></orgPdbQuery>' % (keyword.replace('<', '&lt;').replace('>', '&gt;'))
        response = requests.post('https://www.rcsb.org/pdb/rest/search?sortfield=rank%20Descending', data=data, headers={'Content-Type': 'application/x-www-form-urlencoded'})
        pdbids = response.content.decode().split('\n')
        pdbids = [x for x in pdbids if len(x) > 0]
        self.cached_ids[keyword] = pdbids
        return pdbids

    def count(self, keyword: str) -> int:
        return len(self.query_all(keyword))

    def search(self, keyword: str, rng: range) -> str:
        # Query IDs
        pdbids = self.query_all(keyword)
        pdbids = pdbids[rng.start : rng.stop]

        # Query titles
        form = {
            'pdbids': ','.join(pdbids),
            'customReportColumns': 'structureId,structureTitle',
            'format': 'csv'
        }
        response = requests.post('http://www.rcsb.org/pdb/rest/customReport.csv', form)

        mappings = response.content.decode().replace('<br />', '\n').split('\n')[1:]  # Remain compatible when upstream fixes this bug
        c = csv.reader(mappings)
        results = []
        for (pdbid, title) in zip(pdbids, (x[1] for x in c)):
            result = Result(
                pdbid,
                'https://www.rcsb.org/structure/' + pdbid,
                title
            )
            results.append(result)
        return results

class UniProtSource(Source):
    """
    Source: https://www.uniprot.org/
    """
    pass


class TaxonomySource(Source):
    """
    Source: https://www.uniprot.org/taxonomy/
    """
    pass


class GeneExpressionSource(Source):
    """
    Source: https://www.ncbi.nlm.nih.gov/geo/
    """
    pass


if __name__ == '__main__':
    import sys

    mapping = {
        'human-gene': HumanGeneSource,
        'igem-part': iGEMPartSource,
        'rcsb': RCSBSource,
        'uniprot': UniProtSource,
        'taxonomy': TaxonomySource,
        'gene-expr': GeneExpressionSource
    }

    source_name = sys.argv[1]
    keyword = sys.argv[2]

    for (k, v) in mapping.items():
        if k.startswith(source_name):
            source = v()
            break
    else:
        raise RuntimeError('Unknown data source')

    total = source.count(keyword)
    results = source.search(keyword, range(0, total))
    print('total %d results' % total)
    for (i, result) in enumerate(results):
        print('%d. %s (%s) %s' % (1+i, result.title, result.url, result.summary))
