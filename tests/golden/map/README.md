# Map compatibility fixtures

Anonymous BED inputs and complete outputs were captured with bedtools2.31.1
(`/opt/homebrew/Cellar/bedtools/2.31.1/bin/bedtools`, SHA256
`fc994e1fcfdce8ce8591630d71f21fe67e2275c2a8d967dc0e3dfb4cfa4677a3`).
The original packets retain argv, stdout, stderr, status and input identities:

- `/Volumes/KIOXIA/Developments/tmp/bed-map-contract-20261102.6g3jZT`,
  SEALED.json SHA256 `fc09ee8a9e57ceff27d1d72925567274a4b62f8cab280ecb78709b4f3bf5f1ee`.
- `/Volumes/KIOXIA/Developments/tmp/bed-map-sweep-supplement-20261103.gECDtG`,
  SEALED.json SHA256 `73c73c48149334a76894aae9e6cc5399320aee77ae152ec26268a83c80fb5220`.

The header, nonlexical/lexical, boundary-point and origin-point fixtures come
from the second packet. Other fixtures come from the first. Output bodies here
are unchanged source observations, not regenerated native expectations.

`equal-start-corrected.stdout` is the source invocation with corrected argument
cardinality; it still contains the source's A point-cache defect. The native
test explicitly corrects only point-after to four hits/sum17
(`BEDMAP-2311-POINT-CACHE-LOSS`). Both B equal-start orders must yield the literal
point/count1/sum7 result in the native test
(`BEDMAP-2311-POINT-LOOKAHEAD-LOSS`); the original ordinary-before-point source
observation is retained in the second packet, where it incorrectly has no hit.
These are deliberate repairs, not claims of byte compatibility for those cases.

Numeric-all23 uses column5 and the23 operations in `tests/map.rs`. Text uses
column4. `custom-controls` uses columns4,5 with distinct,mean, delimiter`|`,
null`NA`, precision4. `broadcast-selected` uses columns4,5/count with same strand,
A fraction0.1/B fraction0.25/either. Ordinary split uses collapse,count on
columns4,5. Other option sets are explicit in the complete-output tests.
