"""Regenerate one function at a time, keeping dense suites bounded in memory."""
from collections import defaultdict
from itertools import zip_longest
from functools import lru_cache
from .case_csv import read_cases

from . import ROOT


@lru_cache(maxsize=1)
def profiles():
    result=defaultdict(list)
    for case in read_cases(ROOT/'cases/recipes'):result[case['api']].append(case)
    return dict(result)


class CaseSet:
    def __init__(self,original,count,seed):
        from .explicit_cases import APIS
        self.original=defaultdict(list)
        for case in original:self.original[case['api']].append(case)
        self.count,self.seed=count,seed
        self.apis=sorted(set(self.original)|set(profiles())|set(APIS))

    def function(self,api):
        from .registered_cases import expanded_registered
        from .explicit_cases import cases
        batch=list(self.original.get(api,[]))
        if api in profiles():batch.extend(expanded_registered(self.count,self.seed,apis=[api]))
        batch.extend(cases(self.count,self.seed,apis=[api]))
        return sorted(batch,key=lambda c:c['id'])

    def for_module(self,module):
        for api in self.apis:
            if api.split('.')[0]==module:yield from self.function(api)

    def module_counts(self,module):
        from .explicit_cases import APIS
        names=[a for a in self.apis if a.split('.')[0]==module]
        count=0
        for api in names:
            count+=len(self.original.get(api,[]))
            if api in profiles():count+=len(profiles()[api])+2*self.count
            if api in APIS:count+=2*self.count
        return len(names),count

    def __iter__(self):
        for api in self.apis:yield from self.function(api)

    def __len__(self):
        from .explicit_cases import APIS
        return sum(map(len,self.original.values()))+sum(map(len,profiles().values()))+2*self.count*(len(profiles())+len(APIS))

    def __eq__(self,other):return all(a==b for a,b in zip_longest(self,other))
