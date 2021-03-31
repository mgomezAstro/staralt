### STARALT Example

```python
from astropy.coordinates import SkyCoord
from staralt.object_visibility import ObjVisibility

obj_list = ["PN G226.7+05.6", "PN G326.6+42.2", "PN G011.1+07.0", "PN G208.5+33.2"]

coords = []
for obj in obj_list:
	coords.append(SkyCoord.from_name(obj))
ra = [c.ra.value for c in coords]
dec = [c.dec.value for c in coords]

ov = ObjVisibility(date_obs="2021-4-10", location="spm", ra=ra, dec=dec, names=obj_list, saveFig=True, formatFig="eps")
ov.staralt()
```