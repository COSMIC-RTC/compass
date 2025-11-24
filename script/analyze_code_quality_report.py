import json


def analyze_code_quality_report(file_path):
    # Load the JSON data from the file
    with open(file_path, "r") as file:
        data = json.load(file)

    # Initialize a dictionary to hold the count of issues by severity
    issues_by_severity = {}

    # Iterate through each issue in the JSON data
    for issue in data:
        severity = issue.get("severity", "unknown")
        # Increment the count for the current severity
        issues_by_severity[severity] = issues_by_severity.get(severity, 0) + 1
        location = issue.get("location", "NoLocation")
        path = location.get("path", "NoPath")
        line = location.get("lines", "NoLines")
        if line != "NoLines":
            line = line.get("begin", "NoBegin")
        location = f"{path}:{line}"
        description = issue.get("description", "NoDescription")
        if severity == "blocker":
            print(f"Blocker issue found: {location} {description}")
        elif severity == "critical":
            print(f"Critical issue found: {location} {description}")
        # elif severity == 'major':
        #     print(f"Major issue found: {issue.get('message', 'No message')}")
        # elif severity == 'minor':
        #     print(f"Minor issue found: {issue.get('message', 'No message')}")
        # elif severity == 'info':
        #     print(f"Info issue found: {issue.get('message', 'No message')}")
        # else:
        #     print(f"Unknown severity issue found: {issue.get('message', 'No message')}")

    # Print the summary of issues by severity
    print("Summary of issues by severity:")
    for severity, count in issues_by_severity.items():
        print(f"{severity.capitalize()}: {count}")


# Assuming the JSON data is stored in 'gl-code-quality-report.json'
file_path = "gl-code-quality-report.json"

analyze_code_quality_report(file_path)
