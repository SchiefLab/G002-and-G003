# Getting Data from AWS

We store our data on [AWS S3 Buckets](https://aws.amazon.com/s3/). We have two buckets, one for G002 and one for G003. To get access to the buckets, you need to obtain credentials. Buckets are HIPAA compliant and are encrypted.

# Obtaining AWS credentials.

To get AWS credentials, you can either email [Jordan](mailto:jwillis@scripps.edu) or [Troy](mailto:tsincomb@iavi.org) to get access.

There are two scenarios for access.

1. You already have an AWS key. In that case, email us your IAM ARN and we will add you to the IAM group.

2. You don't have an AWS account and need to be added to the SchiefLab group. In that case, you will receive an email with login instructions, your OTP (one time password) and your security credentials. The security credentials will be your AWS key and secret key and will be used to access the data programatically

## AWS G003

To get G003 data, we use AWS buckets. The bucket `S3URI` is `s3://iavig003sabucket/g003/`. The data is organized by sequencing, sorting and output. To get the data, you need the AWS CLI.

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
$ aws s3 ls s3://iavig003sabucket/ --region af-south-1
PRE g003/
PRE raw_sequences
```

</div>
!!! info "G003"

    If you are shown the `PRE g003/` output, you have access to the bucket.

### AWS S3 Copy Specific Files from AWS

Below is a way to only copy certain components of the AWS bucket

=== "Sorting"

    <div class="termy">

    ```console
    $ aws s3 cp --recursive  s3://iavig003sabucket/g003/sorting/ buckets/g003/sorting
    ---> 100%
    ```
    </div>

=== "Sequencing"

    Get all sequencing files
    <div class="termy">
    ```console
    $ aws s3 cp --recursive  s3://iavig003sabucket/g003/sequencing/ buckets/g003/sequencing
    ---> 100%
    ```
    </div>

    get the sequencing and sorting directory excluding large bcl and fastq files
    <div class="termy">
    ```console
    $ aws s3 cp --recursive s3://iavig003sabucket/g003/ ./g003/sequencing --exclude *working_directory/* --exclude *.fastq.gz --exclude *.tif --exclude *.cbcl --exclude *.imf1 --exclude *.filter --exclude *.bin --exclude *Logs/* --exclude *_stdout --exclude *_stderr --exclude *Autofocus/* --exclude *Intensities/*
    ---> 100%
    ```
    </div>

=== "Output"

    <div class="termy">
    ```console
    $ aws s3 cp --recursive  s3://iavig003sabucket/g003/output/ buckets/g003/output
    ---> 100%
    ```
    </div>

### AWS S3 Copy Local Files to AWS

Copy a directory to AWS. This will copy the entire directory and all subdirectories.

<div class="termy">
```console
$ aws s3 cp --recursive  221118_VH00124_107_AAAW2V3HV  s3://iavig003sabucket/raw_sequnces/
---> 100%
```
</div>

!!! info "221118_VH00124_107_AAAW2V3HV"

    This is the flow cell directory from the sequencer. If you are an end user, update to `s3://iavig003sabucket/raw_sequences`.

### AWS S3 Sync All Files from AWS

SYNC THE ENTIRE BUCKET!!

<div class="termy">
```console
$ aws s3 sync --delete  s3://iavig003sabucket/ buckets/ --region af-south-1
---> 100%
```
</div>

!!! warning "Warning: Large File"

    The entire bucket will likely be over 2 TB

!!! warning "--delete"

    The `--delete` flag will delete any files in the destination that are not in the source. This is useful for keeping the destination in sync with the source. If you don't want to delete files, remove the `--delete` flag.

!!! info "Sync"

    The `sync` command can also be used on specific file susbsets, e.g. sorting

!!! info "s3://iavig003sabucket/g003"

    The subpath `g003` is the main g003 bucket. The top level may have some other app related things.
