from webdav4.client import Client
import os

class NextCloudClient(object):
    """Convenience wrapper class for a webdav4 client used to download files from nextcloud.

    Will take URL, username and password from environment variables:
    NEXTCLOUD_USERNAME
    NEXTCLOUD_PASSWORD and
    NEXTCLOUD_URL
    """

    def __init__(self):

        # Read username and pasword from environment variables
        NEXTCLOUD_USERNAME = os.environ["NEXTCLOUD_USERNAME"]
        NEXTCLOUD_PASSWORD = os.environ["NEXTCLOUD_PASSWORD"]
        NEXTCLOUD_URL = os.environ["NEXTCLOUD_URL"]

        self.client = Client(NEXTCLOUD_URL, auth=(NEXTCLOUD_USERNAME, NEXTCLOUD_PASSWORD))

    def listdir(self, directory):
        """Will list all files in directory on the nextcloud server."""

        return [file_dict['name'] for file_dict in self.client.ls(directory)]

    def download_file(self, file_path, local_file):
        """Will download a file from the nextcloud server and save it as local_file."""
        try:
            self.client.download_file(file_path, local_file)
        except Exception as e:
            print(f"Could not download {file_path} because of: {e}")