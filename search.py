#!/usr/bin/env python3

from typing import List, Dict
import requests
import json
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
        self.jsessionid = {}
        self.qrid = {}
        self.total = {}
        self.cached_results = {}

        # keyword is in direct if and only if the keyword redirects
        # directly to a structure.  This happens when the keyword
        # itself is a PDB ID.
        self.direct = set()

    def get_qrid(self, keyword: str) -> str:
        if keyword in self.qrid:
            return self.qrid[keyword]
        query_string = urlencode({'f': '', 'q': keyword})
        url = 'http://www.rcsb.org/pdb/search/navbarsearch.do?' + query_string
        response = requests.get(url, allow_redirects=False)
        redirect_url = response.headers['Location']
        if redirect_url.startswith('http://www.rcsb.org/pdb/search/structidSearch.do'):
            self.direct.add(keyword)
            return None
        d = parse_qs(response.headers['Location'])
        tmp = d['jsessionid'][0].split('?')
        self.jsessionid[keyword] = tmp[0]
        self.qrid[keyword] = tmp[1].split('=')[1]

        return self.qrid[keyword]

    def count(self, keyword: str) -> int:
        qrid = self.get_qrid(keyword)

        if keyword in self.direct:
            return 1

        if keyword in self.total:
            return self.total[keyword]
        
        qs = urlencode({
            'tabtoshow': 'Current',
            'qrid': qrid,
            'startAt': '0',
            'resultsperpage': '1'
        })
        url = 'http://www.rcsb.org/pdb/json/searchresults.do?' + qs
        response = requests.get(url)
        d = json.loads(response.content.decode())
        self.total[keyword] = d['Pageable']['Total Results Count']
        return self.total[keyword]

    def query_all(self, keyword: str) -> int:
        if keyword in self.cached_results:
            return self.cached_results[keyword]

        if keyword in self.direct:
            url = 'https://www.rcsb.org/structure/' + keyword
            response = requests.get(url)
            soup = BeautifulSoup(response.content, features='html.parser')
            summary = soup.find(id='structureTitle').text
            self.cached_results[keyword] = [
                Result(
                    keyword,
                    url,
                    summary
                )
            ]
            return self.cached_results[keyword]
        
        total = self.count(keyword)
        qrid = self.get_qrid(keyword)
        qs = urlencode({
            'tabtoshow': 'Current',
            'qrid': qrid,
            'startAt': '0',
            'resultsperpage': str(total)
        })
        url = 'http://www.rcsb.org/pdb/json/searchresults.do?' + qs
        response = requests.get(url)
        d = json.loads(response.content.decode())
        if 'Result Set' not in d:
            return []

        results = []
        for t in d['Result Set']:
            result = Result(
                t['PDB ID'],
                'http://www.rcsb.org/structure/' + t['PDB ID'],
                t['Title']
            )
            results.append(result)

        self.cached_results[keyword] = results
        return results

    def search(self, keyword: str, rng: range) -> List[Result]:
        results = self.query_all(keyword)
        return results[rng.start : rng.stop]

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
