#!/usr/bin/env python

import dxpy
import sys
import math
import itertools

P = ""

def IN():
  global P
  P = P + ": "

def OUT():
  global P
  P = P[:-2]

def c(left, right, description):
  if left != right:
    print "ERROR: " + description + " are different:"
    print left
    print right
    sys.exit(1)

def uniq(items):
  single_items = []
  seen = set()
  for i in items:
    if i not in seen:
      seen.add(i)
      single_items.append(i)
  return single_items

def mask(item):
  if type(item) == dict:
    if "$dnanexus_link" in item:
      return "link-xxx"
    else:
      return {key: mask(value) for (key, value) in item.items()}
  elif type(item) == list:
    return [mask(i) for i in item]
  else:
    return item

excluded = frozenset(["project", "id", "created", "modified", "createdBy", "folder", "links", "sponsored"])
social = frozenset(["name", "tags", "properties"])
def compare_descriptions(d1, d2, compare_social, mask_links):
  print P+"Comparing descriptions for objects",d1['id'],"and",d2['id']
  for key in uniq(d1.keys() + d2.keys()):
    if ((key in excluded) or ((not compare_social) and (key in social))):
      continue
    val1 = d1.get(key)
    val2 = d2.get(key)
    if (key == "details") and mask_links:
      val1 = mask(val1)
      val2 = mask(val2)
    c(val1, val2, "describe values for '" + key + "'")
  print P+"Descriptions for objects",d1['id'],"and",d2['id'],"match!"

def compare_objects(o1, o2, compare_social = False, mask_links = True):
  print P+"Comparing object",o1.get_id(),"to",o2.get_id()
  IN()
  d1 = o1.describe(True, True)
  d2 = o2.describe(True, True)

  if d1['id'] == d2['id']:
    OUT()
    print P+"Objects",o1.get_id(),"and",o2.get_id(),"are identical"
    return

  compare_descriptions(d1, d2, compare_social, mask_links)

  # Additional known type comparisons
  if "TrackSpec" in d1['types']:
    compare_trackspecs(o1, o2)
  elif "Wiggle" in d1['types']:
    compare_wiggles(o1, o2)
  elif d1['class'] == "file":
    compare_files(o1, o2)
  elif d1['class'] == "gtable":
    compare_gtables(o1, o2)
  OUT()
  print P+"Objects",o1.get_id(),"and",o2.get_id(),"are identical"

def compare_gtables(gtable1, gtable2):
  d1 = gtable1.describe(True, True)
  d2 = gtable2.describe(True, True)

  g1 = gtable1.iterate_rows(want_dict=True)
  g2 = gtable2.iterate_rows(want_dict=True)

  milestones = set([math.trunc(d1['length']*x/10) for x in range(1,11)])

  print P+"Comparing gtable contents between",d1['id'],"and",d2['id']
  print P,
  i = 0
  for row1 in g1:
    row2 = g2.next()
    for key in row1.keys(): 
        c(row1[key], row2[key], "row contents")

    i = i + 1
    if i in milestones:
      print "%.0f%% " % (i*100.0/d1['length']),
      sys.stdout.flush()
  print
  print P+"Gtables",d1['id'],"and",d2['id'],"have identical contents"

def compare_files(file1, file2):
  print P+"Comparing file contents between",file1.get_id(),"and",file2.get_id()
  while True:
    chunk = file1.read(1024*1024)
    chunk2 = file2.read(1024*1024)
    if len(chunk) == 0 and len(chunk2) == 0:
      break
    if chunk != chunk2:
      c(file1.get_id(), file2.get_id(), "file contents")
  print P+"Files",file1.get_id(),"and",file2.get_id(),"have identical contents"

def compare_trackspecs(o1, o2):
  print P+"Comparing trackspecs between",o1.get_id(),"and",o2.get_id()
  d1 = o1.describe(True, True)
  d2 = o2.describe(True, True)
  for i, (bpp_threshold, rendering_spec) in enumerate(d1['details']['representations']):
    t = rendering_spec['type']
    if t == "empty":
      continue
    s1 = dxpy.get_handler(rendering_spec['source'])
    s2 = dxpy.get_handler(d2['details']['representations'][i][1]['source'])
    
    IN()
    if t == "mappings" or t == "variants" or t == "genes" or t == "spans":
      if (s1.get_id() == o1.get_id()) and (s2.get_id() == o2.get_id()):
        compare_gtables(s1, s2)
      else:
        compare_objects(s1, s2)
    elif t == "wiggle":
      if (s1.get_id() == o1.get_id()) and (s2.get_id() == o2.get_id()):
        compare_wiggles(s1, s2)
      else:
        compare_objects(s1, s2)
    else:
      print "Unknown trackspec type " + t
      sys.exit(1)
    OUT()
  print P+"Trackspecs",o1.get_id(),"and",o2.get_id(),"are identical"

def compare_wiggles(wiggle1, wiggle2):
  print P+"Comparing wiggle",wiggle1.get_id(),"to",wiggle2.get_id()
  d1 = wiggle1.describe(True, True)
  d2 = wiggle2.describe(True, True)

  sources1 = [x['source']['$dnanexus_link'] for x in d1['details']['signals']]
  sources2 = [x['source']['$dnanexus_link'] for x in d2['details']['signals']]

  IN()
  for source in itertools.izip(sources1, sources2):
    compare_gtables(dxpy.get_handler(source[0]), dxpy.get_handler(source[1]))
  OUT()
  print P+"Wiggles",wiggle1.get_id(),"and",wiggle2.get_id(),"are identical"

##############

@dxpy.entry_point('main')
def main(obj1, obj2, compare_social = False, mask_links = True):
  compare_objects(dxpy.get_handler(obj1), dxpy.get_handler(obj2), compare_social, mask_links)
  return {}

dxpy.run()
