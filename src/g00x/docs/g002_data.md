# Contents

- [Contents](#contents)
- [Summary](#summary)
  - [Simple Setup](#simple-setup)
    - [Getting Data from AWS](#getting-data-from-aws)
    - [Obtaining AWS credentials.](#obtaining-aws-credentials)
    - [AWS G002](#aws-g002)
    - [AWS S3 Copy Specific Files](#aws-s3-copy-specific-files)
    - [AWS S3 Sync All Files](#aws-s3-sync-all-files)
  - [Advanced Setup](#advanced-setup)
  - [AWS Instance Setup](#aws-instance-setup)
    - [Box Auto Sync](#box-auto-sync)
    - [Globus Auto Sync](#globus-auto-sync)
  - [Local Machine Setup](#local-machine-setup)
    - [Mounting Box](#mounting-box)
    - [Rclone](#rclone)
    - [Mounting Box](#mounting-box-1)
    - [Mounting Globus](#mounting-globus)

# Summary

AWS s3 buckets contain a stable version of the data. The data is organized by sequencing, sorting and output. To get the data using [Simple Setup](#simple-setup), you need the AWS CLI.

If you are a developer and would like the most recent version of the data with syncing functionality, please see [Advanced Setup](#advanced-setup).

## Simple Setup

### Getting Data from AWS

We store our data on [AWS S3 Buckets](https://aws.amazon.com/s3/). We have two buckets, one for G002 and one for G003. To get access to the buckets, you need to obtain credentials. Buckets are HIPAA compliant and are encrypted.

### Obtaining AWS credentials.

To get AWS credentials, you can either email [Jordan](mailto:jwillis@scripps.edu) or [Troy](mailto:tsincomb@iavi.org) to get access.

There are two scenarios for access.

1. You already have an AWS key. In that case, email us your IAM ARN and we will add you to the IAM group.

2. You don't have an AWS account and need to be added to the SchiefLab group. In that case, you will receive an email with login instructions, your OTP (one time password) and your security credentials. The security credentials will be your AWS key and secret key and will be used to access the data programatically

### AWS G002

To get G002 data, we use AWS buckets. The bucket `S3URI` is `s3://iavig002westbucket/g002/`. The data is organized by sequencing, sorting and output. To get the data, you need the AWS CLI.

To get AWS CLI, follow the instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

Once you get the AWS CLI, you need to configure it. Follow the instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html).

<div class="termy">

```console
$ aws configure
AWS Access Key ID [None]: Secret key in email credentials
AWS Secret Access Key [None]: Secret key in email credentials
Default region name [None]: us-west-2
Default output format [None]: json
```

</div>

Once you are configured, we can check if you have access to the S3 bucket.

<div class="termy">

```console
$ aws s3 ls s3://iavig002westbucket/
PRE g002/
```

</div>
!!! info "G002"

    If you are shown the `PRE g002/` output, you have access to the bucket.

### AWS S3 Copy Specific Files

Below is a way to only copy certain components of the AWS bucket

=== "Sorting"

    <div class="termy">

    ```console
    $ aws s3 cp --recursive  s3://iavig002westbucket/g002/G002/sorting/ buckets/g002/G002/sorting
    ---> 100%
    ```
    </div>

=== "Sequencing"

    Get all sequencing files
    <div class="termy">
    ```console
    $ aws s3 cp --recursive  s3://iavig002westbucket/g002/G002/sorting/ buckets/g002/G002/sequencing
    ---> 100%
    ```
    </div>

    get the sequencing directory excluding large bcl and fastq files
    <div class="termy">
    ```console
    $ aws s3 cp --recursive s3://iavig002westbucket/g002/G002/sequencing ./g002/G002/sequencing --exclude *working_directory/* --exclude *.fastq.gz --exclude *.tif --exclude *.cbcl --exclude *.imf1 --exclude *.filter --exclude *.bin --exclude *Logs/* --exclude *_stdout --exclude *_stderr --exclude *Autofocus/* --exclude *Intensities/*
    ---> 100%
    ```
    </div>

=== "Output"

    <div class="termy">
    ```console
    $ aws s3 cp --recursive  s3://iavig002westbucket/g002/G002/output/ buckets/g002/G002/output
    ---> 100%
    ```
    </div>

### AWS S3 Sync All Files

<div class="termy">
```console
$ aws s3 sync --delete  s3://iavig002westbucket/ buckets/
---> 100%
```
</div>

!!! warning "Warning: Large File"

    The entire bucket will likely be over 2 TB

!!! warning "--delete"

    The `--delete` flag will delete any files in the destination that are not in the source. This is useful for keeping the destination in sync with the source. If you don't want to delete files, remove the `--delete` flag.

!!! info "Sync"

    The `sync` command can also be used on specific file susbsets, e.g. sorting

## Advanced Setup

For those users that need to mount [Box](https://www.box.com/) or [Globus](https://www.globus.org/), we a have 2 sets of instructions below. The first is for an quick setup if you are using an AWS instance; otherwise for a local machine, we recommend using the second set of instructions.

- The 2 Available Options:
  1. [AWS Instance Quick Setup](#aws-instance-setup)
  2. [Local Machine Setup](#local-machine-setup)


## AWS Instance Setup

For this to work you must have all the following completed:

1. You have a running AWS instance with a home dictory to a user with sudo privileges
2. You were invited to the G002 [Box](https://www.box.com/).
3. You were invited to the G002 [Globus](https://www.globus.org/) collection.
4. G00x is installed on your AWS instance. If not, follow the instructions [here](index.md)

### Box Auto Sync

This could take 3+ hours, but it will be running in the background using systemd.

<div class="termy">

```console
$ g00x g002 box setup
## click the link and sign into your box account when prompted
```

</div>

To check on the status of the sync, run the following below.


<div class="termy">

```console
$ g00x g002 box status
```

</div>

### Globus Auto Sync

This could take 1+ hours, but it will be running in the background using systemd.

<div class="termy">

```console
$ g00x g002 globus setup
## click the link and sign into your globus account when prompted
```

</div>

To check on the status of the sync, run the following below.


<div class="termy">

```console
$ g00x g002 globus status
```

</div>

!!! INFO "You are done! Please skip the rest of this page regarding local machine setup."

## Local Machine Setup

!!! warning "DO NOT DO IF YOU ALREADY COMPLETE THE AWS QUICK SETUP"

### Mounting Box

### Rclone

To mount box, we will use a utility called [R clone](https://rclone.org/). R clone is a command line utility that can be used to mount cloud storage. To install R clone, follow the instructions [here](https://rclone.org/install/).

Once you install, you can run the following

<div class="termy">

```console
$ rclone config
```

</div>

It will ask you a few questions. I used the following answers

<div class="termy">

```console
n) new remote
name > box
storage > 6
client_id > (leave blank)
client_secret > (leave blank)
box_config_file > (leave blank)
box_sub_type > 1
Edit advanced config? (y/n) > n
Use auto config? (y/n) > <see below>
```

</div>

`use auto config` depends on if you are using a machine with a browser (like your personal computer). If you are on a server without access to a browser, hit `n` and you will be given a link to copy and paste into a browser.

If you hit `n`, head over to a machine with rclone installed (e.g. brew install rclone) that has a browser and type

<div class="termy">

```console
$ rclone authorize box
```

</div>

That will take you into a login page. Just login and hit authorize. Come back to the command line console and you will see a code that looks like the following.

```
#Paste this code into the rclone configuration:
>result  {"access_token":"LGOqxaTWf2Tzc6Na0","token_type":"bearer","refresh_token":"t8DBslxjlZ7JXiRBk5rqv1cVFN6O5BUcFEzzzFS2","expiry":"2022-12-19T12:32:17.317533-06:00"}
```

There you now have Box setup to mount! Check `~/.config/rclone/rclone.conf` to see the configuration and that the token is in there.

### Mounting Box

To mount box, you can use the following command

<div class="termy">

```console
you can't mount this on fsx so put on /mnt/
$ mkdir /mnt/box

change ownership
$ sudo -R chown jwillis:jwillis /mnt/box

mount with rclone
$ rclone mount --daemon g002: /mnt/box
```

</div>

Now box will be mounted to /mnt/box

!!! warning "Path must exist"

    The /mnt/path/to/box must exist. If it doesn't, you will get an error.

### Mounting Globus

To install globus, we need Globus Personal Connect. You can find the complete installation [here](https://docs.globus.org/how-to/globus-connect-personal-linux/)

Once you install, you can run the following

<div class="termy">

```console
$ /path/to/globus/install/globusconnectpersonal -start &
```

</div>
