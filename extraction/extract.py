"""Execute one read-only raw-source query using standard PostgreSQL settings.

Connection parameters are resolved by libpq (PGSERVICE / PGHOST / PGPORT /
PGDATABASE / PGUSER and a locally configured password file). No connection
string or patient values are printed. CSV outputs belong in private storage.
"""
import argparse
import csv
from pathlib import Path
import re


def load_queries(dataset):
    text = (Path(__file__).parent / 'sql' / (dataset + '.sql')).read_text('utf-8')
    parts = re.split(r'^-- name: ([a-z0-9_]+)\s*$', text, flags=re.M)
    return {parts[i]: parts[i+1].strip() for i in range(1, len(parts), 2)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset', choices=['mimiciv', 'eicu'], required=True)
    p.add_argument('--query', help='Named query; use --list to inspect names')
    p.add_argument('--list', action='store_true')
    p.add_argument('--stay-ids', type=Path, help='Private CSV with icu_stay_id column')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    queries = load_queries(args.dataset)
    if args.list:
        print('\n'.join(queries))
        return
    if args.query not in queries or args.output is None:
        p.error('Supply a valid --query and --output')
    query = queries[args.query]
    params = None
    if '%(stay_ids)s' in query:
        if args.stay_ids is None:
            p.error('This query requires --stay-ids')
        with args.stay_ids.open(encoding='utf-8-sig', newline='') as f:
            ids = sorted({int(row['icu_stay_id']) for row in csv.DictReader(f)})
        if not ids:
            p.error('No ICU stay identifiers supplied')
        params = {'stay_ids': ids}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    import psycopg2
    with psycopg2.connect('') as conn:
        conn.set_session(readonly=True)
        with conn.cursor() as settings:
            settings.execute("SET statement_timeout TO '30min'")
        with conn.cursor(name='core_extraction') as cur:
            cur.itersize = 10000
            cur.execute(query, params)
            rows = cur.fetchmany(10000)
            # Named PostgreSQL cursors expose column metadata after FETCH.
            with args.output.open('x', encoding='utf-8', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([c.name for c in cur.description])
                total = 0
                while rows:
                    writer.writerows(rows)
                    total += len(rows)
                    rows = cur.fetchmany(10000)
    print(f'Exported {total} source rows. No patient values were logged.')


if __name__ == '__main__':
    main()
