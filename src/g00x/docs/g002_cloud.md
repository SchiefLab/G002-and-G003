# Cloud management

At this moment we are using terraform to describe the cloud instances we are requisitioning.

## Terraform

Terraform is a tool for building, changing, and versioning infrastructure safely and efficiently. Terraform can manage existing and popular service providers as well as custom in-house solutions.

### Installation

For full installation instructions, go to [https://learn.hashicorp.com/tutorials/terraform/install-cli](https://learn.hashicorp.com/tutorials/terraform/install-cli)


#### Mac

For mac use homebrew to install

<div class="termy">
```bash
$ brew tap hashicorp/tap
$ brew install hashicorp/tap/terraform
$ brew update
```
</div>

#### Starting an instance

To start an instance, you need to have a valid AWS account and have the credentials in your environment. You can use the following command to start an instance.

<div class="termy">
```bash
$cd terraform/oregon
$ terraform init
$ terraform apply
```
</div>

!!! note

    This will start an instance as recorded in single_ec2.tf. That is a structured terraform file that describes the resources you need:

    ```

    # tells terraform you are using AWS
    terraform {
    required_providers {
        aws = {
        source  = "hashicorp/aws"
        version = "~> 3.27"
        }
    }

    required_version = ">= 0.14.9"
    }

    # Tells which profile to use in your credentials and to start the instance in Oregon
    provider "aws" {
    profile = "default"
    region  = "us-west-2" # oregon
    }


    # associate a static IP address with your instnace
    resource "aws_eip_association" "eip_assoc" {
    instance_id         = aws_instance.g00x_dedicated_instance.id # can assocaite this with an instnace or network interface
    public_ip           = "44.240.169.113"  #elastic IP
    allow_reassociation = true
    }



    # This is the instance you want to start
    resource "aws_instance" "g00x_dedicated_instance" {
        ami = "ami-0fd14da38e402236e" # Ubuntu LTS 20.04 - with Luster Kernel and CellRanger

        # what size of the instnace
        instance_type = "m5.12xlarge"
        key_name      = "G00x"
        network_interface {
            network_interface_id  = "eni-06717892d1883626c" #network interface, logical grouping of vpcid, subnet and security for G00x
        device_index          = 0
        delete_on_termination = false
        }
        root_block_device {
        delete_on_termination = true
        volume_size           = 500
        tags                  = { Name = "G00x Root Volume" }
        }

        # this is the startup file we will use
        user_data = file("startup.sh") # file directive can install stuff
        tags = {
        Name = "G00x Oregon"
        }
        }

    resource "aws_eip" "example" {
    vpc = true
    }
    ```


#### Stopping an instance

To stop a terraformed instance, we can simply destroy it. This will not delete anything since we are using an distributed FSx file system.

<div class="termy">
```bash
$ terraform destroy
```
</div>
