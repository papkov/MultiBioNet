import requests
import sys


def get_ensembl_info(ensg):
    ext = "/archive/id/%s?" % ensg
    r = requests.get("http://rest.ensembl.org%s" % ext,
                     headers={"Content-Type": "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    return repr(decoded)
