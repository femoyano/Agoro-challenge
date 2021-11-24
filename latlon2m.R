# Get m from lat lon

latlon2m <- function(lat1, lon1, lat2, lon2) {  # generally used geo measurement function
  R = 6378.137 #Radius of earth in KM
  dLat = lat2 * pi / 180 - lat1 * pi / 180
  dLon = lon2 * pi / 180 - lon1 * pi / 180
  a = sin(dLat/2) * sin(dLat/2) +
    cos(lat1 * pi / 180) * cos(lat2 * pi / 180) *
    sin(dLon/2) * sin(dLon/2)
  c = 2 * atan2(sqrt(a), sqrt(1-a))
  d = R * c
  return (d * 1000)# meters
}

latlon2m(43, -90, 40,-97)

# Calc m - somethigns not ok here
mlat <- function(lat) {111320 * lat}
mlon <- function(lat,lon) {lon * 40075000 * cos(lat) / 360}
lat1 = 40
lat2 = 43
lon1 = -90
lon2 = -97
mlat(lat1)
mlat(lat2)
mlon(lat1, lon1)
mlon(lat1, lon2)
