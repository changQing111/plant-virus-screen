"""Versioned, conservative RNA-virus discovery policy. Scores are heuristics."""
from collections import defaultdict

POLICY_VERSION = '1.0.0'
# Stable Pfam accessions; explicitly exclude cellular RDR and auxiliary domains.
# This initial catalytic-core panel is deliberately incomplete, not a pan-virus set.
VIRAL_RDRP = {
    'PF00680': 'RdRP_1', 'PF00978': 'RdRP_2',
    'PF00998': 'RdRP_3', 'PF02123': 'RdRP_4',
    'PF07925': 'RdRP_5', 'PF05919': 'Mitovir_RNA_pol',
}

def allowed_model(accession, name):
    return VIRAL_RDRP.get(accession.split('.')[0]) == name

def regional_competition(rows, overlap=0.5, margin=0.10):
    """Compare local HSPs on the same strand; never add overlapping scores.

    Require reciprocal overlap >= 50%; scores within 10% are unresolved.
    Non-overlapping viral/cellular support is mixed, not cellular exclusion.
    """
    by_query = defaultdict(list)
    details = []
    for r in rows:
        by_query[r['query_id']].append(r)
    summaries = {}
    for query, hits in by_query.items():
        states = set()
        for h in hits:
            cls = h['classification']
            if cls not in ('viral', 'cellular'):
                states.add('unclassified')
                continue
            start, end = int(h['qstart']), int(h['qend'])
            lo, hi = sorted((start, end))
            peers = []
            for p in hits:
                if p['classification'] == cls:
                    continue
                ps, pe = int(p['qstart']), int(p['qend'])
                pl, ph = sorted((ps, pe))
                shared = max(0, min(hi, ph) - max(lo, pl) + 1)
                if (end >= start) == (pe >= ps) and shared / max(hi-lo+1, ph-pl+1) >= overlap:
                    peers.append(p)
            competitor = max(peers, key=lambda p: float(p['bitscore'])) if peers else None
            score = float(h['bitscore'])
            other = float(competitor['bitscore']) if competitor else 0.0
            delta = (score-other)/max(score, other, 1.0)
            state = cls if competitor is None or delta > margin else 'ambiguous' if abs(delta) <= margin else 'outscored'
            if competitor and competitor['classification'] == 'unclassified' and state != cls:
                state = 'ambiguous'
            if state != 'outscored':
                states.add(state)
            details.append(dict(query_id=query, subject_id=h['subject_id'], qstart=start,
                                qend=end, classification=cls, bitscore=score,
                                competing_subject=competitor['subject_id'] if competitor else '',
                                competing_bitscore=other, relative_margin=delta, region_status=state))
        if {'viral','cellular'} <= states:
            status = 'mixed'
        elif 'ambiguous' in states:
            status = 'ambiguous'
        elif 'unclassified' in states:
            status = 'unclassified'
        elif 'viral' in states:
            status = 'viral'
        elif 'cellular' in states:
            status = 'cellular'
        else:
            status = 'unclassified'
        summaries[query] = status
    return summaries, details

def discovery_priority(evidence, competition, mapping_ok, curated=False):
    kinds = set(evidence.split(';'))
    if competition in ('mixed','ambiguous','cellular'):
        return 'review_conflict', 'regional_conflict_requires_review'
    if curated and 'RdRP_domain' in kinds:
        return ('novelty_review_high' if mapping_ok else 'novelty_review_pending_mapping',
                'curated_RdRp_without_nt' if 'nt' not in kinds else 'curated_RdRp_with_nt')
    return 'homology_review', 'homology_or_legacy_domain_evidence'
