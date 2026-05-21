import csv
import urllib.request
import re

ecs = [
    "1.16.3.1", "1.4.7.1", "1.5.1.7", "1.8.1.8", "2.1.1.13", "2.1.1.148",
    "2.1.1.37", "2.3.1.129", "2.4.2.45", "2.4.2.9", "2.6.1.1", "2.6.1.11",
    "2.7.1.50", "2.7.4.7", "2.7.6.3", "2.7.9.1", "3.1.3.15", "3.2.1.22",
    "3.2.1.23", "3.2.1.4", "3.2.1.55", "3.4.11.5", "3.4.24.55", "3.6.1.13",
    "3.6.3.14", "4.1.1.37", "4.1.1.41", "4.2.1.8", "4.3.3.6", "5.1.3.3",
    "5.1.3.n2", "5.4.99.62", "6.1.1.5", "6.1.1.9", "6.3.1.5", "6.3.4.15"
]

results = []

for ec in ecs:
    url = f"http://rest.kegg.jp/get/ec:{ec}"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            content = response.read().decode('utf-8')
            
            name = ""
            ec_class = ""
            
            for line in content.split('\n'):
                if line.startswith('NAME'):
                    # The KEGG NAME field can have multiple lines, or just take the first
                    name_match = re.match(r'^NAME\s+(.+)$', line)
                    if name_match:
                        name = name_match.group(1).strip()
                elif line.startswith('CLASS'):
                    class_match = re.match(r'^CLASS\s+(.+)$', line)
                    if class_match:
                        ec_class = class_match.group(1).strip()
            
            print(f"[{ec}] {name} | {ec_class}")
            results.append({"EC": ec, "Name": name, "Class": ec_class})
    except Exception as e:
        print(f"Error fetching {ec}: {e}")
        results.append({"EC": ec, "Name": "ERROR", "Class": "ERROR"})

with open("results/tables/section_3_3/ec_kegg_annotations.csv", 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=["EC", "Name", "Class"])
    writer.writeheader()
    writer.writerows(results)
