#!/usr/bin/env python

import pyliftover

lo = LiftOver('hg17', 'hg18')
lo = LiftOver('hg17ToHg18.over.chain.gz')

pyliftover.LiftOver()

# FROM: https://github.com/konstantint/pyliftover/tree/master/pyliftover
# convert_coordinate(self, chromosome, position, strand='+'):
#         '''
#         Returns a *list* of possible conversions for a given chromosome position.
#         The list may be empty (no conversion), have a single element (unique conversion), or several elements (position mapped to several chains).
#         The list contains tuples (target_chromosome, target_position, target_strand, conversion_chain_score),
#         where conversion_chain_score is the "alignment score" field specified at the chain used to perform conversion. If there
#         are several possible conversions, they are sorted by decreasing conversion_chain_score.

