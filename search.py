#!/usr/bin/env python3

from typing import List


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
    def search(self, keyword: str) -> List[Result]:
        """
        Search a certain data source for a specific keyword.

        :return: A list of [Result]s.
        """
        raise NotImplementedError()


class HumanGeneSource(Source):
    """
    Source: https://ghr.nlm.nih.gov/gene
    """
    pass


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

    for result in source.search(keyword):
        print('%s: %s' % (result.title, result.url))

