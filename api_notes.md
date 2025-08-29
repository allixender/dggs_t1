
## Comparison of similar work

### current rhealpix functions

We initialize a DGGS instance with an ellipsoid definition, and decide if we work in degrees or radians. Cell/Zone IDs (SUID) are awkward tuples of base face character (N, P, O ..) and children 1-9 int numbers. Also all operations offer to work on plane => False or True

I mostly work with: rddgs = rhealpixdggs.dggs.WGS84_003 and following methods:

- cell from suid: rdggs.cell(suid)
- cell from point: rdggs.cell_from_point()
- get outline:
	cell.boundary(n=2,plane=False)
	cell.vertices(plane=False)
- get centroid:
	cell.centroid()

- region filler:
	cells = rdggs.cells_from_region(res, se, nw, plane=False) (returns a nested 	list indicating topology)

- get parent or children is implicit adding or removing int 1-9, but there is also a programmatical way via centroid and recreating cell at lower resolution
	- coords ... cell.centroid(plane=False)
	  rdggs.cell_from_point(coords..., resolution, plane=False)
	- cell.subcells()

- cell.neighbor, cell.neighbors

Surround Australia additional "glue":

- https://github.com/surroundaustralia/rhealpix-geo/blob/master/rheal/dggs_classes.py
- https://github.com/surroundaustralia/rhealpix-sf/blob/master/rhealsf/dggs_functions.py


### Eaggr API

https://raw.githubusercontent.com/riskaware-ltd/open-eaggr/master/Documents/Programmer's%20Guide.pdf

EAGGR officially support ISEA3H and ISEA4T DGGS, but in practice it seems only ISEA4T is fully functional. Eaggr uses an accuracy concept, which doesn't seem to allow direct access to a resolution/level configuration

- you need to initialse the DGGS system instance
	dggs = Eaggr(Model.ISEA4T)
	dggs.convert_point_to_dggs_cell(LatLongPoint(latitude, longitude. accuracy)

- converts a single point to a DGGS cell: convert_point_to_dggs_cell()
- DggsCell object to a LatLongPoint: convert_dggs_cell_to_point()
- DGGS domain to points in the latitude/longitude domain: convert_dggs_cells_to_points()
- exports the geometry of the DGGS cell as a WKT or GeoJSON string as a polygon: convert_dggs_cell_outline_to_shape_string()
- returns the parent cell(s) for a supplied cell: get_dggs_cell_parents()
- returns the child cells for a supplied cell: get_dggs_cell_children()
- returns the sibling cells for a supplied cell in the same resolution, that share a parent: get_cell_siblings()
- DggsCells and returns the smallest cell that contains all of the supplied cells: get_bounding_dggs_cell()


### H3 base API

H3 is fixed, fixed coded Faces, indices and rotation

- geoToH3(lat, lng, res) ⇒ H3Index
- h3ToGeo(h3Index) ⇒ Array.<number>
- h3ToGeoBoundary(h3Index, [formatAsGeoJson]) ⇒ Array.<Array.<number>>

- kRing(h3Index, ringSize) ⇒ Array.<H3Index>
- kRingDistances(h3Index, ringSize) ⇒ Array.<Array.<H3Index>>
- hexRing(h3Index, ringSize) ⇒ Array.<H3Index>
- polyfill(coordinates, res, [isGeoJson]) ⇒ Array.<H3Index>

- h3GetResolution(h3Index) ⇒ number
- h3ToParent(h3Index, res) ⇒ H3Index
- h3ToChildren(h3Index, res) ⇒ Array.<H3Index>
- h3ToCenterChild(h3Index, res) ⇒ H3Index

- getRes0Indexes() ⇒ Array.<H3Index>
- h3IsResClassIII(h3Index) ⇒ boolean
- getPentagonIndexes(res) ⇒ Array.<H3Index>
- h3IsPentagon(h3Index) ⇒ boolean
- h3GetBaseCell(h3Index) ⇒ number
- h3GetFaces(h3Index) ⇒ Array.<number>
- cellArea(h3Index, unit) ⇒ number
- exactEdgeLength(edge, unit) ⇒ number
- hexArea(res, unit) ⇒ number
- edgeLength(res, unit) ⇒ number
- numHexagons(res) ⇒ number

- h3IndexesAreNeighbors(origin, destination) ⇒ boolean
- compact(h3Set) ⇒ Array.<H3Index>
- uncompact(compactedSet, res) ⇒ Array.<H3Index>
- h3Distance(origin, destination) ⇒ number


### Google S2 base API

S2 is quite advanced and can do a lot geometry stuff independent of DGGS (but its mathematics are on a sphere not ellipsoid)
- it works out of the box with a reliable default orientation
- it seems to be possible to retrieve base faces
- it seems to be possible to initialise a custom setup and rotate/flip the origin/cube of DGGS

- Convert from Lat / Lng:  latLngToKey(lat, lng, level)
- Convert between Hilbert Curve Quadtree Key and S2 Cell Id: keyToId(key)
- and reverse: idToKey(id)
- example: base 4 quadkey to base 10 (decimal) S2 Cell Id: '4/032212303102210' becomes '9749618446378729472'
- Convert between Quadkey and Id: keyToLatLng(key) and idToLatLng(id)
- Neighbours: latLngToNeighborKeys(lat, lng, level)

Python S2:

- S2LatLng.FromDegrees(51.3368602, 0.4931979)
	london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    cell = s2.S2CellId(london)
	polygon = s2.S2Polygon(cell)
	.ToPoint()
	.parent(10)

- cell = s2.S2Cell(s2.S2CellId.FromToken(cell_id)
	cell.GetS2LatLngVertex(pos 0 - 4)
- rect = s2.S2LatLngRect(ll, ur)
- S2RegionCoverer().GetCovering(rect)
- S2 points (not cells yet) can be in-memory indexed: S2PointIndex


## DGGRID

although not meant as an API, main functionality revolves around:

- region filler
	- grid_cell_polygons_for_extent(): fill extent/subset with cells at given 	resolution (clip or world)
	- grid_cellids_for_extent(): get_all_indexes/cell_ids for dggs at given resolution (clip or world)
- generate cells from either IDs or points
	- grid_cell_polygons_from_cellids(): geometry_from_cellid for dggs at given 	resolution (from id list): cell_ids to cells
	- cells_for_geo_points(): poly_outline for point/centre at given resolution

- sample/aggregate point data into dggs cells, i.e. binning
