#!/usr/bin/env python3

from typing import List, Tuple
import requests
from bs4 import BeautifulSoup
from urllib.parse import urlencode


class Result:
    def __init__(self, title, url, summary):
        self.title = title
        self.url = url
        self.summary = summary


class Source:
    """
    An abstract base class for all data sources.

    The only required method is `search`.
    """
    def search(self, keyword: str, offset: int) -> Tuple[int, List[Result]]:
        """
        Search a certain data source for a specific keyword.

        :return: (total number of results, list of [Result]s)
        """
        raise NotImplementedError()


class HumanGeneSource(Source):
    """
    Source: https://ghr.nlm.nih.gov/gene
    """
    def search(self, keyword: str, offset: int = 0):
        params = {
            'query': keyword,
            'start': offset
        }
        query_string = urlencode(params)
        response = requests.get('https://ghr.nlm.nih.gov/search?%s' % query_string)
        soup = BeautifulSoup(response.content, features='html.parser')
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
        total = int(soup.find(class_='ss-item').text[5:-1])
        return (total, results)


class iGEMPartSource(Source):
    """
    Source: http://parts.igem.org/Main_Page
    """
    pass


class RCSBSource(Source):
    """
    Source: http://www.rcsb.org/
    """
    pass


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

    total, results = source.search(keyword)
    print('total %d results' % total)
    for result in results:
        print('%s: %s' % (result.title, result.url))

