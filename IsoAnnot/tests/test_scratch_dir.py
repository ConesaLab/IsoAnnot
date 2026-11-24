import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from scratch_dir import filesystem_type, is_network_filesystem, make_scratch_dir

# /proc/mounts of a cluster node: /home and /storage on Ceph, a local disk for scratch
MOUNTS = ("/dev/sda2 / ext4 rw,relatime 0 0\n"
          "tmpfs /tmp tmpfs rw,nosuid,nodev 0 0\n"
          "/dev/nvme0n1 /scr_local xfs rw,noatime 0 0\n"
          "10.2.0.5,10.2.0.6:/storage /storage ceph rw,relatime,name=hpc 0 0\n"
          "10.2.0.5,10.2.0.6:/storage/home /home ceph rw,relatime,name=hpc 0 0\n"
          "server:/export /mnt/my\\040data nfs4 rw,relatime 0 0\n"
          "ceph-fuse /cephfuse fuse.ceph rw,nosuid 0 0\n")


@pytest.fixture
def mounts_file(tmp_path):
    mounts = tmp_path / "mounts"
    mounts.write_text(MOUNTS)
    return str(mounts)


@pytest.mark.parametrize("path, fs_type", [
    ("/scr_local/patienza", "xfs"),
    ("/home/patienza/isoannot_trial/out", "ceph"),
    ("/storage/gge/Pablo", "ceph"),
    ("/storagex", "ext4"),          # not inside /storage: only the / mount contains it
    ("/tmp", "tmpfs"),
    ("/mnt/my data/run", "nfs4"),
    ("/cephfuse/run", "fuse.ceph"),
    ("/usr/local", "ext4"),
])
def test_filesystem_type(mounts_file, path, fs_type):
    """The type is the one of the longest mount point that contains the path."""
    assert filesystem_type(path, mounts_file) == fs_type


def test_filesystem_type_without_mounts_file(tmp_path):
    """Without /proc/mounts (not Linux) the type is unknown, and unknown is not taken as a network filesystem."""
    assert filesystem_type("/tmp", str(tmp_path / "missing")) is None
    assert not is_network_filesystem(None)


@pytest.mark.parametrize("fs_type, network", [
    ("ceph", True), ("fuse.ceph", True), ("nfs", True), ("nfs4", True), ("lustre", True), ("gpfs", True),
    ("ext4", False), ("xfs", False), ("tmpfs", False), ("overlay", False), ("fuseblk", False),
])
def test_is_network_filesystem(fs_type, network):
    assert is_network_filesystem(fs_type) == network


def test_make_scratch_dir(tmp_path, mounts_file):
    """A new folder for each job inside scratch_dir, which is created if it is missing."""
    scratch = tmp_path / "scratch" / "patienza"
    first = make_scratch_dir(str(scratch), "isoannot_repeatmasker", mounts_file)
    second = make_scratch_dir(str(scratch), "isoannot_repeatmasker", mounts_file)
    assert first != second
    for work_dir in (first, second):
        assert os.path.isdir(work_dir)
        assert os.path.dirname(work_dir) == str(scratch)
        assert os.path.basename(work_dir).startswith("isoannot_repeatmasker_")


def test_make_scratch_dir_default(tmp_path, mounts_file, monkeypatch):
    """Without scratch_dir the folder is made in the system's temporary folder."""
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    monkeypatch.setattr("tempfile.tempdir", None)
    work_dir = make_scratch_dir("", "isoannot_repeatmasker", mounts_file)
    assert os.path.dirname(work_dir) == str(tmp_path)


def test_make_scratch_dir_stops_on_network_filesystem(tmp_path, caplog):
    """A scratch folder on Ceph stops the run before RepeatMasker starts, and nothing is created in it."""
    scratch = tmp_path / "home" / "patienza" / "scratch"
    scratch.mkdir(parents=True)
    mounts = tmp_path / "mounts"
    mounts.write_text(f"/dev/sda2 / ext4 rw 0 0\n10.2.0.5:/storage/home {tmp_path}/home ceph rw 0 0\n")
    with pytest.raises(SystemExit):
        make_scratch_dir(str(scratch), "isoannot_repeatmasker", str(mounts))
    assert f"The scratch folder {scratch} is on a network filesystem (ceph)." in caplog.text
    assert "--scratch" in caplog.text
    assert os.listdir(scratch) == []
