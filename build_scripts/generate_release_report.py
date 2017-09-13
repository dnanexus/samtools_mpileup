# Get commits in issues
from __future__ import unicode_literals, print_function
import requests
import argparse
import jinja2
import json
import re
import pytz
from collections import namedtuple
from tempfile import TemporaryFile
from base64 import b64encode
from os import getcwd, remove
from os.path import join, isfile
from datetime import datetime, date
from jira import JIRA
from subprocess import check_call, check_output
from FixedOffset import FixedOffset


GitHubCommit = namedtuple('GitHubCommit', 'origin_user origin_user_link source_repo hash link message is_merged')
JenkinsTest = namedtuple('JenkinsTest', 'recent_build_status recent_build_link recent_build_time recent_build_number issue_obj')
ChangeEntry = namedtuple('ChangeEntry', 'Field Date From To By')

JIRA_ISSUE_RE = re.compile(r'^([A-Z][A-Z]+\-\d+).*')  # string must start with issue key


def create_issue_jenkins_links(jira, jenkins_tests):
    def remote_link(jenkins_test):
        """TODO: When access is granted
        issue=jira.issue('SCI-452')  # Subtask: Check for broken links
        destination = {
            "url":"https://scijenkins.internal.dnanexus.com/view/App(let)%20Test/job/sniffles/482/",
            "title":"SCI-452-job-run",
            "summary":"A job run on this date",
            "icon": {
                "url16x16":"https://scijenkins.internal.dnanexus.com/static/a4ce876b/images/headshot.png",
                "title":"Support Ticket"
                }
        }
        globalId="system={}".format(destination['url'])

        remote_link_obj = jira.add_remote_link(issue, destination, globalId=globalId)
        """
        pass
    pass


def jenkins_recent_test(jenkins_url, user, auth_token, job_name):
    """Gets results of most recent Test execution for each issue. If results are present returns JenkinsTest object"""
    auth = "--user {user}:{token}".format(user=user, token=auth_token)
    curl_cmd = 'curl {server_url}/job/{job_name}/lastBuild/api/json {auth}'.format(
        server_url=jenkins_url, job_name=job_name, auth=auth)
    response = json.loads(check_output(curl_cmd, shell=True))  # Assuming empty dict is returned, need to check
    if not response:
        return
    # GET console output
    build_numb = response['number']
    time_utc_no_offset = datetime.fromtimestamp(int(response['timestamp'] / 1000)).strftime('%Y-%m-%d-(%H:%M:%S)')
    log_name = "{job_name}-{utc_time}".format(job_name=job_name, utc_time=time_utc_no_offset)
    curl_console_cmd = 'curl -o {log_name} "https://scijenkins.internal.dnanexus.com/job/{job_name}/{build_numb}/logText/progressiveText?start=0" {auth}'.format(
        log_name=log_name, job_name=job_name, build_numb=build_numb, auth=auth)
    check_call(curl_console_cmd)
    return JenkinsTest(
        recent_build_status=response['result'] == 'SUCCESS', recent_build_link=response['url'],
        recent_build_time=response['timestamp'], recent_build_number=build_numb, recent_build_console="")


def jenkins_request_maker(jenkins_url, user, auth_token):
    """Return function to make Jenkins API calls"""
    def make_jenkins_request(resource, method="GET"):
        return check_output(curl_cmd.format(
            server_url=jenkins_url, auth=auth, method=method, resource=resource.lstrip('/')), shell=True)
    auth = "--user {user}:{token}".format(user=user, token=auth_token)
    curl_cmd = 'curl -X {method} {server_url}/{resource} {auth}'
    return make_jenkins_request


def instantiate_jira(jira_url, user, pswrd):
    return JIRA(jira_url, basic_auth=(user, pswrd))


def internal_jira_api_requester(jira_url, internal_api_path, user, pswrd):
    """Returns a function to make specified internal JIRA api calls"""
    def make_request(resource, query=None, data=None, req_type='GET'):
        # Construct URL
        url = "{api}{resource}{query}".format(
            api=api_url.rstrip('/'),
            resource="/" + resource,
            query="?" + query if query else "")
        # Construct requests func call
        req_type = req_type.upper()
        req = {
            "POST": requests.post,
            "GET": requests.get,
            "PUT": requests.put
        }.get(req_type)
        input_dict = {"url": url, 'headers': auth_header}
        if req is None:
            raise Exception("{} is not a supported request type".format(req_type))
        if req_type == "POST" and data is not None:
            # TODO logic to determine data or json field in request.post
            input_dict["json"] = data
        # log request url
        resp = req(**input_dict)
        try:
            return resp.json()
        except ValueError:
            return resp

        return resp.json

    auth_header = {"Authorization": "Basic {encoded_cred}".format(
        encoded_cred=b64encode("{user}:{pswrd}".format(user=user, pswrd=pswrd)))
    }
    api_url = join(jira_url, internal_api_path)
    # consider changing make_request.__name__ before return
    return make_request


def _make_html_report(output_html_path, release_name, issues, commits, jenkins_tests):
    """Uses template in current directory"""
    templateLoader = jinja2.FileSystemLoader(searchpath="./")
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "quality_report_template.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    with open(output_html_path, 'w') as out_file:
        out_file.write(template.render(release_name=release_name, issues=issues, commits=commits, jenkins_tests=jenkins_tests))


def subset_changelog(issue, start_dt, end_dt):
    """subset issue changelog to specified dates

    FIXME: DO NOT USE unambiguos datetime objects, get pytz and use that
    INSPIRATION: https://community.atlassian.com/t5/Answers-Developer-Questions/Is-it-possible-to-get-the-issue-history-using-the-REST-API/qaq-p/510094"""
    def datetime_tz_aware(date_str):
        """https://stackoverflow.com/questions/1101508/how-to-parse-dates-with-0400-timezone-string-in-python"""
        naive_date_str, _, offset_str = date_str.rpartition('-')
        offset_str = "-" + offset_str
        naive_dt = datetime.strptime(naive_date_str, '%Y-%m-%dT%H:%M:%S.%f')
        offset = int(offset_str[-4:-2]) * 60 + int(offset_str[-2:])
        if offset_str[0] == "-":
            offset = -offset
        dt = naive_dt.replace(tzinfo=FixedOffset(offset))
        return dt

    changelog = issue.changelog
    subset_logs = []
    # need to convert history.created : datetime.strptime('1984-06-02T19:05:00.000Z', '%Y-%m-%dT%H:%M:%S.%fZ')
    for history in changelog.histories:
        history_dt = datetime_tz_aware(history.created).replace(tzinfo=None)
        if history_dt < start_dt or history_dt > end_dt:  # If the entries are sorted I can break when history_dt > end_dt, check this
            continue
        human_history = history_dt.strftime('%Y-%m-%d-%H:%M')
        for item in history.items:
            try:

                subset_logs.append(ChangeEntry(
                    Field=item.field, Date=human_history,
                    From=item.fromString if item.fromString else "None",
                    To=item.toString if item.toString else "None",
                    By=history.author.displayName))
            except TypeError:
                print('Field failed: ', item.field)

    return subset_logs


def generate_release_report(jira, dev_status_req, jenkins_make_req, release_id, overwrite=False):
    """Generate release report for provided release id.
    TODO: Add various templates for reports and expose selection
    TODO: Double check internal request will always have detail[0],commit, & repo fields
    """
    def _parse_resp(response):
        repo_objs = response['detail'][0]['repositories']  # Double check no KeyError possible with this request
        commits = []
        for repo in repo_objs:
            for commit in repo['commits']:
                commits.append(GitHubCommit(
                    origin_user=commit['author']['name'], origin_user_link=commit['author']['url'],
                    source_repo=repo['name'], hash=commit['id'], link=commit['url'],
                    message=commit['message'], is_merged=commit['merge']))
        return commits

    def _issue_match(job_name, issues):
        """Find first matching job."""
        match = JIRA_ISSUE_RE.match(job_name)
        if not match:
            return False
        target_key = match.group(1)
        matching_issue = next((issue for issue in issues if issue.key == target_key), False)
        if matching_issue:
            return {'matching_issue': matching_issue, 'job_name': job_name}

    def create_jenkins_test_obj(issue_job_pair):
        """Get build and test results from test page
        Add remotelink to test from issue: https://developer.atlassian.com/jiradev/jira-platform/guides/other/guide-jira-remote-issue-links
        """
        response = json.loads(jenkins_make_req(resource="job/{job_name}/lastCompletedBuild/api/json".format(job_name=issue_job_pair['job_name'])))
        return JenkinsTest(
            recent_build_status=response['result'] == 'SUCCESS', recent_build_link=response['url'],
            recent_build_time=datetime.fromtimestamp(int(response['timestamp'] / 1000)).strftime('%Y-%m-%d-(%H:%M:%S)'),
            recent_build_number=response['number'], issue_obj=issue_job_pair['matching_issue'])

    def _version_start_date():
        try:  # time in 'YYYY-MM-DD'
            start = datetime.strptime(version.startDate, "%Y-%m-%d").replace(tzinfo=None)
        except AttributeError:
            print('WARNING: Specify a Release version start date in JIRA')
            raise Exception('Release version has no official start')

        try:
            end = datetime.strptime(version.releaseDate, "%Y-%m-%d").replace(tzinfo=None)
        except AttributeError:
            end = datetime.utcnow()

        return start, end

    report_path = join(getcwd(), 'Quality-Report-{}.html'.format(date.today().isoformat()))
    if isfile(report_path):
        if overwrite:
            remove(report_path)
        else:
            print("Report already exist generation skipped")
            return

    version = jira.version(release_id)
    start_utc_datetime, end_utc_datetime = _version_start_date()
    project = jira.project(version.projectId)
    jql_str = "project = {proj} AND fixVersion = {release_name}".format(
        proj=project.name, release_name=version.name)
    release_issues = jira.search_issues(jql_str, expand='changelog')
    for issue in release_issues:
        response = dev_status_req(resource='issue/detail', query="issueId={}&applicationType=github&dataType=repository".format(issue.id))
        issue.commits = _parse_resp(response)
        issue.subsethistory = subset_changelog(issue, start_dt=start_utc_datetime, end_dt=end_utc_datetime)

    job_list = [job['name'] for job in json.loads(jenkins_make_req("api/json"))['jobs']]  # arr job names
    print('JOB LIST: ', job_list)
    existing_issue_jobs = [_issue_match(job, release_issues) for job in job_list]
    existing_issue_jobs = [item for item in existing_issue_jobs if item]  # {'matching_issue': matching_issue, 'job_name': job_name}

    jenkins_tests = [create_jenkins_test_obj(issue_job) for issue_job in existing_issue_jobs]  # Consider uploading and attaching logs to issue here
    _make_html_report(
        output_html_path=report_path,
        release_name=version.name,
        issues=release_issues,
        commits=[commit for issue in release_issues for commit in issue.commits],
        jenkins_tests=jenkins_tests)

    # create_issue_jenkins_links(jira_requester, j_tests)  # Create remote link to jenkins


def _parse_args():
    parser = argparse.ArgumentParser(description='Generate Release Version Report',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('jira_url', help="Base URL for jira installation")
    parser.add_argument('release_id', help="Release ID, if provided generates a version release report.")

    parser.add_argument('-u', '--user', help="JIRA username", default="sosazuwa@dnanexus.com", type=str, dest='user')
    parser.add_argument('-p', '--pass', help="JIRA password", type=str, dest='pswrd')
    parser.add_argument('--jenkins-url', help="Jenkins url", type=str, dest='jenkins_url')
    parser.add_argument('--jenkins-auth', help="Jenkins Auth token", type=str, dest='jenkins_auth')
    parser.add_argument("--overwrite", help="Overwrites old files with generated files", action="store_true", dest="overwrite")
    return parser.parse_args()


def main(args):
    if not args.release_id:
        # TODO validate release ID
        print("Please specify a valid release ID")
        exit(1)  # Use correct non-zero exit
    generate_release_report(
        dev_status_req=internal_jira_api_requester(
            jira_url=args.jira_url, internal_api_path='rest/dev-status/1.0', user=args.user, pswrd=args.pswrd),
        jenkins_make_req=jenkins_request_maker(args.jenkins_url, user=args.user, auth_token=args.jenkins_auth),
        jira=instantiate_jira(args.jira_url, user=args.user, pswrd=args.pswrd),
        release_id=args.release_id, overwrite=args.overwrite)


if __name__ == '__main__':
    args = _parse_args()
    main(args)
