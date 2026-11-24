#!/usr/bin/env python3
"""
Makes the temporary folder a rule works in and prints its path. Used by rule run_repeatmasker.

RepeatMasker writes thousands of small temporary files. On a network filesystem (Ceph, NFS, Lustre...) this
can block the filesystem for every user of a cluster, so the folder has to be on a local disk: the run
stops here when it is not, before RepeatMasker starts.
"""

import argparse, logging, os, sys, tempfile

# Filesystem types shared over the network, as /proc/mounts names them ("fuse.ceph" is read as "ceph")
NETWORK_FILESYSTEMS = {"afs", "beegfs", "ceph", "cifs", "glusterfs", "gpfs", "lustre", "nfs", "nfs4", "panfs",
                       "pvfs2", "smb3", "sshfs"}


def filesystem_type(path, mounts_file="/proc/mounts"):
    """
    Type of the filesystem path is on: the one of the longest mount point that contains it.

    :return: the type as written in mounts_file, or None when there is no mounts_file (not Linux).
    """
    if not os.path.isfile(mounts_file):
        return None
    path = os.path.realpath(path)
    best_mount, best_type = "", None
    with open(mounts_file) as mounts:
        for line in mounts:
            fields = line.split()
            if len(fields) < 3:
                continue
            # Spaces in a mount point are written as \040
            mount_point, fs_type = fields[1].replace("\\040", " "), fields[2]
            inside = path == mount_point or path.startswith(mount_point.rstrip("/") + "/")
            if inside and len(mount_point) >= len(best_mount):
                best_mount, best_type = mount_point, fs_type
    return best_type


def is_network_filesystem(fs_type):
    return fs_type is not None and fs_type.removeprefix("fuse.") in NETWORK_FILESYSTEMS


def make_scratch_dir(scratch_dir, prefix, mounts_file="/proc/mounts"):
    """
    Create a new folder for one job inside scratch_dir, or inside the system's temporary folder when
    scratch_dir is empty. The caller removes it.

    :return: path of the new folder. Stops the run if scratch_dir is on a network filesystem.
    """
    scratch_dir = scratch_dir or tempfile.gettempdir()
    os.makedirs(scratch_dir, exist_ok=True)
    fs_type = filesystem_type(scratch_dir, mounts_file)
    if is_network_filesystem(fs_type):
        logging.error(f"The scratch folder {scratch_dir} is on a network filesystem ({fs_type}). RepeatMasker writes "
                      "thousands of temporary files, which can block a shared filesystem for every user. "
                      "Give a folder on a local disk with --scratch (on a cluster, the local disk of the "
                      "compute nodes), or leave the layer out with --config repeat_masker=False. Abort!")
        sys.exit(1)
    return tempfile.mkdtemp(prefix=f"{prefix}_", dir=scratch_dir)


def main():
    parser = argparse.ArgumentParser(description="Make a temporary work folder on a local disk and print its path")
    parser.add_argument("--scratch_dir", default="", help="Folder to make it in (default: the system's temporary folder)")
    parser.add_argument("--prefix", default="isoannot", help="Start of the name of the new folder (default: %(default)s)")
    args = parser.parse_args()
    print(make_scratch_dir(args.scratch_dir, args.prefix))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
