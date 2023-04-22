def clean_text(regex, text):
    import re
    found_matches = re.sub(regex, "\n", text)
    return found_matches

def search_text(regex, text):
    import re
    rgx = re.compile(regex)
    found_matches = re.findall(rgx, text)
    return found_matches

def main():
    import glob, os
    files = [y for x in os.walk('.') for y in glob.glob(os.path.join(x[0], '*.m'))]
    for f in files:
        with open(f, 'r') as content_file:
            content = content_file.read()
        content = clean_text(" +\n", content)
        with open(f, 'w') as content_file:
            content_file.write(content)

if __name__ == "__main__":
    main()

