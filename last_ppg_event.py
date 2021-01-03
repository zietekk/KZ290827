import arcpy
from math import sqrt, acos, pi
import numpy as np
import pandas as pd

def create_point_feature_class(pnts):
    arcpy.CreateFeatureclass_management("C:/Users/ziete/Desktop/ppg_egz/wyniki", "Concave_Hull.shp", "POINT")
    cur = arcpy.da.InsertCursor("C:/Users/ziete/Desktop/ppg_egz/wyniki/Concave_Hull.shp", ["SHAPE@"])

    points = []
    for i in pnts:
        for j in i:
            points.append(j)

    for row in points:
        cur.insertRow([row])

def concave_hull(input_feature_class, output_feature_class, k_0=3, field_choice="", includenull=True):
    try:
        import arcpy
        import itertools
        import math
        import os
        import sys
        import traceback
        import string

        arcpy.overwriteOutput = True

        # Funkcja zwraca liste OID dla najblizszych k
        def kNeighbours(k, oid, pDict, excludeList=[]):
            hypotList = [math.hypot(pDict[oid][0] - pDict[id][0], pDict[oid][1] - pDict[id][1]) for id in pDict.keys() if id <> oid and id not in excludeList]
            hypotList.sort()
            hypotList = hypotList[0:k]
            oidList = [id for id in pDict.keys() if math.hypot(pDict[oid][0] - pDict[id][0], pDict[oid][1] - pDict[id][1]) in hypotList and id <> oid and id not in excludeList]
            return oidList

        # Funkcja zwraca liste X,Y po obroceniu pktu wokol innego pktu
        def RotateXY(x, y, xc=0, yc=0, angle=0):
            x = x - xc
            y = y - yc
            xr = (x * math.cos(angle)) - (y * math.sin(angle)) + xc
            yr = (x * math.sin(angle)) + (y * math.cos(angle)) + yc
            return [xr, yr]

        # Funcja znajduje OID elementu dla prawego kata
        def Rightmost(oid, angle, pDict, oidList):
            origxyList = [pDict[id] for id in pDict.keys() if id in oidList]
            rotxyList = []
            for p in range(len(origxyList)):
                rotxyList.append(RotateXY(origxyList[p][0], origxyList[p][1], pDict[oid][0], pDict[oid][1], angle))
            minATAN = min([math.atan2((xy[1] - pDict[oid][1]), (xy[0] - pDict[oid][0])) for xy in rotxyList])
            rightmostIndex = rotxyList.index([xy for xy in rotxyList if math.atan2((xy[1] - pDict[oid][1]), (xy[0] - pDict[oid][0])) == minATAN][0])
            return oidList[rightmostIndex]

        # Funkcja wykrywania samoczynnego przeciecia sie warstw
        def selfIntersects(polyline):
            lList = []
            selfIntersects = False
            for n in range(0, len(line.getPart(0)) - 1):
                lList.append(arcpy.Polyline(arcpy.Array([line.getPart(0)[n], line.getPart(0)[n + 1]])))
            for pair in itertools.product(lList, repeat = 2):
                if pair[0].crosses(pair[1]):
                    selfIntersects = True
                    break
            return selfIntersects

        # Funkcja budujaca
        def createHull(pDict, outCaseField, lastValue, kStart, dictCount, includeNull):
            # k obejmuje wszystkie pkty danych
            enclosesPoints = False
            notNullGeometry = False
            k = kStart

            if dictCount > 1:
                pList = [arcpy.Point(xy[0], xy[1]) for xy in pDict.values()]
                mPoint = arcpy.Multipoint(arcpy.Array(pList), sR)
                minY = min([xy[1] for xy in pDict.values()])

                while not enclosesPoints and k <= 30:
                    # Znajduje poczatkowy pkt o najmniejszej wartosci wsp Y
                    startOID = [id for id in pDict.keys() if pDict[id][1] == minY][0]
                    # Znajduje nastepny pkt w prawo
                    kOIDList = kNeighbours(k, startOID, pDict, [])
                    minATAN = min([math.atan2(pDict[id][1] - pDict[startOID][1], pDict[id][0] - pDict[startOID][0]) for id in kOIDList])
                    nextOID = [id for id in kOIDList if math.atan2(pDict[id][1] - pDict[startOID][1], pDict[id][0] - pDict[startOID][0]) == minATAN][0]
                    # Tablica granic
                    bArray = arcpy.Array(arcpy.Point(pDict[startOID][0], pDict[startOID][1]))
                    bArray.add(arcpy.Point(pDict[nextOID][0], pDict[nextOID][1]))

                    currentOID = nextOID
                    prevOID = startOID

                    excludeList = [startOID, nextOID]

                    steps = 2
                    while currentOID <> startOID and len(pDict) <> len(excludeList):
                        try:
                            angle = math.atan2((pDict[currentOID][1] - pDict[prevOID][1]), (pDict[currentOID][0] - pDict[prevOID][0]))
                            oidList = kNeighbours(k, currentOID, pDict, excludeList)
                            nextOID = Rightmost(currentOID, 0 - angle, pDict, oidList)
                            pcArray = arcpy.Array([arcpy.Point(pDict[currentOID][0], pDict[currentOID][1]), arcpy.Point(pDict[nextOID][0], pDict[nextOID][1])])
                            while arcpy.Polyline(bArray, sR).crosses(arcpy.Polyline(pcArray, sR)) and len(oidList) > 0:
                                excludeList.append(nextOID)
                                oidList.remove(nextOID)
                                oidList = kNeighbours(k, currentOID, pDict, excludeList)
                                if len(oidList) > 0:
                                    nextOID = Rightmost(currentOID, 0 - angle, pDict, oidList)
                                    pcArray = arcpy.Array([arcpy.Point(pDict[currentOID][0], pDict[currentOID][1]),
                                                           arcpy.Point(pDict[nextOID][0], pDict[nextOID][1])])
                            bArray.add(arcpy.Point(pDict[nextOID][0], pDict[nextOID][1]))
                            prevOID = currentOID
                            currentOID = nextOID
                            excludeList.append(currentOID)
                            steps += 1
                            if steps == 4:
                                excludeList.remove(startOID)
                        except ValueError:
                            arcpy.AddMessage("Nie znaleziono najblizszych sasiadow z " + str(pDict[currentOID]) + " , rozszerzenie poszukiwania")
                            break
                    # Zamkniecie granicy, sprawdzenie
                    bArray.add(arcpy.Point(pDict[startOID][0], pDict[startOID][1]))
                    pPoly = arcpy.Polygon(bArray, sR)
                    if pPoly.length == 0:
                        break
                    else:
                        notNullGeometry = True
                    if mPoint.within(arcpy.Polygon(bArray, sR)):
                        enclosesPoints = True
                    else:
                        k += 1
                if not mPoint.within(arcpy.Polygon(bArray, sR)):
                    arcpy.AddWarning("Punkty odstajace, dane nie sa zamkniete")

            # Wstawienie poligonow
            if (notNullGeometry and includeNull == False) or includeNull:
                if outCaseField > " ":
                    insFields = [outCaseField, "POINT_CNT", "ENCLOSED", "SHAPE@"]
                else:
                    insFields = ["POINT_CNT", "ENCLOSED", "SHAPE@"]
                rows = arcpy.da.InsertCursor(outFC, insFields)
                row = []
                if outCaseField > " ":
                    row.append(lastValue)
                row.append(dictCount)
                if notNullGeometry:
                    row.append(enclosesPoints)
                    row.append(arcpy.Polygon(bArray, sR))
                else:
                    row.append(-1)
                    row.append(None)
                rows.insertRow(row)
                del row
                del rows
            elif outCaseField > " ":
                arcpy.AddMessage("\nGeometria nie moze byc pusta dla przypadku " + str(lastValue))
            else:
                arcpy.AddMessage("\nGeometria nie moze byc pusta")

        # Klasa wejsciowa
        inPoints = input_feature_class
        inDesc = arcpy.Describe(inPoints)
        inPath = os.path.dirname(inDesc.CatalogPath)
        sR = inDesc.spatialReference


        k = int(k_0)
        kStart = k

        # Klasa wynikowa
        outFC = output_feature_class
        outPath = os.path.dirname(outFC)
        outName = os.path.basename(outFC)

        caseField = field_choice
        if caseField > " ":
            fields = inDesc.fields
            for field in fields:
                if field.name == caseField:
                    caseFieldType = field.type
                    if caseFieldType not in ["SmallInteger", "Integer", "Single", "Double", "String", "Date"]:
                        arcpy.AddMessage("\nNieprawidlowy typ pola " + caseField + "\n")
                        caseField = " "
                    else:
                        if caseFieldType in ["SmallInteger", "Integer", "Single", "Double"]:
                            caseFieldLength = 0
                            caseFieldScale = field.scale
                            caseFieldPrecision = field.precision
                        elif caseFieldType == "String":
                            caseFieldLength = field.length
                            caseFieldScale = 0
                            caseFieldPrecision = 0
                        else:
                            caseFieldLength = 0
                            caseFieldScale = 0
                            caseFieldPrecision = 0

        # Definicja nazwy pliku wynikowego
        outCaseField = str.upper(str(caseField))
        if outCaseField == "ENCLOSED":
            outCaseField = "ENCLOSED1"
        if outCaseField == "POINT_CNT":
            outCaseField = "POINT_CNT1"
        if outFC.split(".")[-1] in ("shp", "dbf"):
            outCaseField = outCaseField[0:10]

        includeNull = includenull

        inDesc = arcpy.Describe(inPoints)
        sR = inDesc.spatialReference
        arcpy.env.OutputCoordinateSystem = sR
        oidName = str(inDesc.OIDFieldName)
        if inDesc.dataType == "FeatureClass":
            inPoints = arcpy.MakeFeatureLayer_management(inPoints)

        # Tworzenie wynikowej klasy
        if '.SHP' in outName.upper():
            outName = outName[:-4]
        arcpy.AddMessage(outPath + "; " + outName)
        outFC = arcpy.CreateFeatureclass_management(outPath, outName, "POLYGON", "#", "#", "#", sR).getOutput(0)
        if caseField > " ":
            if caseFieldType in ["SmallInteger", "Integer", "Single", "Double"]:
                arcpy.AddField_management(outFC, outCaseField, caseFieldType, str(caseFieldScale),
                                          str(caseFieldPrecision))
            elif caseFieldType == "String":
                arcpy.AddField_management(outFC, outCaseField, caseFieldType, "", "", str(caseFieldLength))
            else:
                arcpy.AddField_management(outFC, outCaseField, caseFieldType)
        arcpy.AddField_management(outFC, "POINT_CNT", "Long")
        arcpy.AddField_management(outFC, "ENCLOSED", "SmallInteger")

        # Budowanie struktur
        rowCount = 0
        caseCount = 0
        dictCount = 0
        pDict = {}  # Slownik z OID, X i Y bez duplikowanych pktow
        if caseField > " ":
            fields = [caseField, 'OID@', 'SHAPE@X', 'SHAPE@Y']
            valueDict = {}
            with arcpy.da.SearchCursor(inPoints, fields) as searchRows:
                for searchRow in searchRows:
                    keyValue = searchRow[0]
                    if not keyValue in valueDict:
                        valueDict[keyValue] = [[searchRow[1], searchRow[2], searchRow[3]]]
                    else:
                        valueDict[keyValue].append([searchRow[1], searchRow[2], searchRow[3]])
            for lastValue in sorted(valueDict):
                caseCount += 1
                for p in valueDict[lastValue]:
                    rowCount += 1
                    if [p[1], p[2]] not in pDict.values():
                        pDict[p[0]] = [p[1], p[2]]
                        dictCount += 1
                createHull(pDict, outCaseField, lastValue, kStart, dictCount, includeNull)
                pDict = {}
                dictCount = 0
        else:
            fields = ['OID@', 'SHAPE@X', 'SHAPE@Y']
            for p in arcpy.da.SearchCursor(inPoints, fields):
                rowCount += 1
                if [p[1], p[2]] not in pDict.values():
                    pDict[p[0]] = [p[1], p[2]]
                    dictCount += 1
                    lastValue = 0
            createHull(pDict, outCaseField, lastValue, kStart, dictCount, includeNull)

        if caseField == " " and arcpy.GetParameterAsText(3) > " ":
            arcpy.AddMessage("\n" + arcpy.GetParameterAsText(3) + " zostalo pominiete ze wzledu na bledna nazwe!")



    # Errory
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)

        msgs = "GP ERRORS:\n" + arcpy.GetMessages(2) + "\n"
        arcpy.AddError(msgs)

    arcpy.CheckGeometry_management(output_feature_class, "C:/Users/ziete/Desktop/ppg_egz/wyniki/spr")

    if arcpy.GetCount_management("C:/Users/ziete/Desktop/ppg_egz/wyniki/spr")[0] == "2":
        return "Error"
    else:
        return polygon_to_polyline(output_feature_class)

def part_geometry(part, feature, check, minimal_geometry_list):
    point_number = -1
    results = []
    for i in part[0:-1]:
        point_number += 1
        point_before = part[0:-1][point_number - 1]
        point = part[0:-1][point_number]
        if point_number != len(part[0:-1]) - 1:
            point_next = part[0:-1][point_number + 1]
        else:
            point_next = part[0:-1][0]

        length_in = segment_length(point_before, point)
        length_out = segment_length(point_next, point)
        angle_in = vertex_angle(point_before, point, point_next, feature, check)

        minimal_distances = []
        for i in minimal_geometry_list:
            if i == "Error":
                minimal_distances.append("Unknown")
            else:
                minimal_distances.append(deflection(point, i))

        results.append([point_number, length_in, length_out, angle_in, minimal_distances[0], minimal_distances[1], minimal_distances[2], minimal_distances[3], minimal_distances[4], minimal_distances[5]])

    return results

def feature_geometry(parts, feature, minimalne_geometrie):
    iterator = 0
    intermediate_results = []
    for i in parts:
        iterator += 1
        if iterator == 1:
            results = part_geometry(i, feature, 'normal', minimalne_geometrie)
            intermediate_results.append(results)
        else:
            anArray = arcpy.Array()
            for j in i:
                anArray.add(j)
            new_feature = arcpy.Polygon(anArray)
            results = part_geometry(i, new_feature, 'reverse', minimalne_geometrie)
            intermediate_results.append(results)

    results_2 = []
    for i in intermediate_results:
        for j in i:
            results_2.append(j)

    return results_2

def polygon_to_polyline(polygon):
    arcpy.PolygonToLine_management(polygon, polygon[0:-4] + "_linia.shp")

    polylines = []
    for row in arcpy.da.SearchCursor(polygon[0:-4] + "_linia.shp", ["SHAPE@"]):
        polylines.append(row[0])

    return polylines[0]

# DLUGOSC
def segment_length(p1, p2):
    x1 = p1.X
    y1 = p1.Y

    x2 = p2.X
    y2 = p2.Y

    dx = x2 - x1
    dy = y2 - y1

    return sqrt(dx**2 + dy**2)

# KAT
def vertex_angle(p1, p2, p3, feature, check):
    a = segment_length(p1, p2)
    b = segment_length(p2, p3)
    c = segment_length(p1, p3)

    if round(a+b, 10) == round(c, 10):
        angle = 180
    else:
        angle = acos((a**2 + b**2 - c**2)/(2*a*b))*180/pi

        x1 = p1.X
        y1 = p1.Y

        x3 = p3.X
        y3 = p3.Y

        middle_point = arcpy.Point()
        middle_point.X, middle_point.Y = (x1 + x3)/2, (y1 + y3)/2

        if check == 'normal':
            if feature.disjoint(middle_point):
                angle = 360 - angle
            else:
                pass
        else:
            if feature.disjoint(middle_point):
                pass
            else:
                angle = 360 - angle

    return angle

# CZY WIERZCHOLEK JEST WEZLEM
def is_node(vertex, building_id, features):
    expression = "gmlID <> " + "'" + building_id + "'"

    count_of_neighbor_vertex = 0
    for row in arcpy.da.SearchCursor(features, ["OID@", "SHAPE@", "gmlId"], where_clause = expression):
        for part in row[1]:
            polygon = arcpy.Polygon(part)
            if intersect(vertex, polygon):
                count_of_neighbor_vertex = count_of_neighbor_vertex + 1

    return count_of_neighbor_vertex

# STRZALKA
def deflection(vertex, polyline):
    return polyline.distanceTo(vertex)

# PRZECIECIE
def intersect(geom1, geom2):
    if not geom1.disjoint(geom2):
        return True
    else:
        return False

# FUNKCJA MINIMALNYCH OCZEK
def minimal_geometry(feature):
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Users/ziete/Desktop/ppg_egz/wyniki/rectangle_by_area.shp", "RECTANGLE_BY_AREA")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Users/ziete/Desktop/ppg_egz/wyniki/rectangle_by_width.shp", "RECTANGLE_BY_WIDTH")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Users/ziete/Desktop/ppg_egz/wyniki/convex_hull.shp", "CONVEX_HULL")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Users/ziete/Desktop/ppg_egz/wyniki/circle.shp", "CIRCLE")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Users/ziete/Desktop/ppg_egz/wyniki/envelope.shp", "ENVELOPE")

    list = [r"C:/Users/ziete/Desktop/ppg_egz/wyniki/rectangle_by_area.shp", r"C:/Users/ziete/Desktop/ppg_egz/wyniki/rectangle_by_width.shp", r"C:/Users/ziete/Desktop/ppg_egz/wyniki/convex_hull.shp", r"C:/Users/ziete/Desktop/ppg_egz/wyniki/circle.shp", r"C:/Users/ziete/Desktop/ppg_egz/wyniki/envelope.shp"]

    polygon_geometry_list = []
    polyline_geometry_list = []
    for i in list:
        polyline_geometry_list.append(polygon_to_polyline(i))
        for row in arcpy.da.SearchCursor(i, ["SHAPE@"]):
            polygon_geometry_list.append(row[0])

    return polyline_geometry_list


def building(input_feature_class, building_id):
    expression = "gmlID = " + "'" + building_id + "'"

    for row in arcpy.da.SearchCursor(input_feature_class, ["OID@", "SHAPE@", "gmlId"], where_clause=expression):
        minimal_geometry_list = minimal_geometry(row[1])
        parts_list = []
        for part in row[1]:
            single_part = []
            for pnt in part:
                if pnt:
                    single_part.append(pnt)
                else:
                    trigger = True
                    parts_list.append(single_part)
                    single_part = []
            if len(parts_list) == 0 or trigger:
                trigger = False
                parts_list.append(single_part)

            create_point_feature_class(parts_list)
            minimal_geometry_list.append(concave_hull("C:/Users/ziete/Desktop/ppg_egz/wyniki/Concave_Hull.shp", "C:/Users/ziete/Desktop/ppg_egz/wyniki/Concave_Hull_wynik.shp"))
            results = feature_geometry(parts_list, arcpy.Polygon(part), minimal_geometry_list)

    for i in results:
        i.insert(0, building_id)

    return results


def main():
    arcpy.AddMessage("\nTrwa przetwarzanie... ")
    arcpy.env.overwriteOutput = 1

    input_feature_class = r"C:/Users/ziete/Desktop/ppg_egz/dane/BUBD.shp"

    building_ids = []
    for row in arcpy.da.SearchCursor(input_feature_class, ["OID@", "SHAPE@", "gmlId"]):
        building_ids.append(row[2])

    intermediate_results = []
    for i in building_ids:
        intermediate_results.append(building(input_feature_class, i))

    results = []
    for i in intermediate_results:
        for j in i:
            results.append(j)

    df = pd.DataFrame(np.array(results),
                      columns = ['id_budynku', 'wierzcholek', 'dlugosc_segm_przed', 'dlugosc_segm_po', 'kat_wewn', 'strzalka_RECTANGLE_BY_AREA', 'strzalka_RECTANGLE_BY_WIDTH', 'strzalka_CONVEX_HULL', 'strzalka_CIRCLE', 'strzalka_ENVELOPE', 'strzalka_CONCAVE_HULL'])

    df.to_csv('C:/Users/ziete/Desktop/results.csv', index = False)
    arcpy.AddMessage("\nZakonczono przetwarzanie")

if __name__ == '__main__':
    main()