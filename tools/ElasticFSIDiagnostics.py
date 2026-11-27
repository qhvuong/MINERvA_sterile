"""ROOT-independent accounting for the elastic-FSI selected-event study."""
import bisect
import heapq
import math
from collections import defaultdict


GROUPS = ('negative', 'zero', 'suppressed', 'near_one', 'enhanced_to_2', 'above_2')


def weight_group(weight, tolerance=1e-8):
    if not math.isfinite(weight):
        raise ValueError('Nonfinite elastic FSI weight')
    if weight < 0:
        return 'negative'
    if weight == 0:
        return 'zero'
    if abs(weight - 1) <= tolerance:
        return 'near_one'
    if weight < 1:
        return 'suppressed'
    return 'enhanced_to_2' if weight <= 2 else 'above_2'


class Accounting:
    def __init__(self, edges, top_n=10):
        self.edges = list(map(float, edges))
        if len(self.edges) < 2 or any(not math.isfinite(x) for x in self.edges):
            raise ValueError('Invalid EN4 edges')
        if any(a >= b for a, b in zip(self.edges, self.edges[1:])) or top_n < 1:
            raise ValueError('Edges must increase and top_n must be positive')
        self.top_n = top_n
        self.sums = defaultdict(lambda: [0, 0., 0., 0.])
        self.tops = defaultdict(list)
        self.zeros = defaultdict(list)
        self.sequence = 0

    def add(self, category, energy, baseline, elastic, context):
        corrected = baseline * elastic
        delta = corrected - baseline
        if not all(math.isfinite(x) for x in (energy, baseline, elastic, corrected, delta)):
            raise ValueError('Nonfinite event: {}'.format(context))
        group = weight_group(elastic)
        # ROOT convention: last upper edge belongs to overflow.
        bin_number = bisect.bisect_right(self.edges, energy)
        region = ('underflow' if bin_number == 0 else
                  'overflow' if bin_number == len(self.edges) else 'regular')
        for label in ('TOTAL', category):
            for key in ((label, region, group), (label, 'bin_{}'.format(bin_number), group)):
                row = self.sums[key]
                row[0] += 1
                row[1] += baseline
                row[2] += corrected
                row[3] += delta
        event = dict(context, category=category, EN4_GeV=energy, bin=bin_number,
                     region=region, w430=baseline, elastic=elastic, w431=corrected,
                     delta=delta, weight_group=group)
        # Rank regular bins separately so overflow cannot obscure the reported shift.
        for ranking in ('absolute_delta', 'elastic'):
            score = abs(delta) if ranking == 'absolute_delta' else elastic
            heap = self.tops[(category, region, ranking)]
            item = (score, self.sequence, event)
            if len(heap) < self.top_n:
                heapq.heappush(heap, item)
            elif item[:2] > heap[0][:2]:
                heapq.heapreplace(heap, item)
        if elastic == 0 and len(self.zeros[(category, region)]) < 2:
            self.zeros[(category, region)].append(event)
        self.sequence += 1

    def rows(self):
        rows = []
        for (category, scope, group), (count, baseline, corrected, delta) in sorted(self.sums.items()):
            rows.append(dict(category=category, scope=scope, weight_group=group,
                             events=count, sum_w430=baseline, sum_w431=corrected,
                             delta=delta, ratio=corrected / baseline if baseline else None))
        return rows

    def totals(self):
        totals = defaultdict(lambda: [0, 0., 0., 0.])
        for (category, scope, group), values in self.sums.items():
            if scope not in ('regular', 'underflow', 'overflow'):
                continue
            for i, x in enumerate(values):
                totals[(category, scope)][i] += x
        return [dict(category=c, scope=s, events=v[0], sum_w430=v[1], sum_w431=v[2],
                     delta=v[3], ratio=v[2]/v[1] if v[1] else None)
                for (c, s), v in sorted(totals.items())]

    def top_events(self):
        return [dict(category=c, scope=s, ranking=r,
                     events=[x[2] for x in sorted(heap, reverse=True)])
                for (c, s, r), heap in sorted(self.tops.items())]
