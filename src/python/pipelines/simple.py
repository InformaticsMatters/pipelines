import sys, argparse

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='Simple example')
    parser.add_argument('--name', help='your name')

    args = parser.parse_args()
    execute(args)

'''
Main execution entry point
'''
def execute(args):
    result = 'Hello ' + args.name
    print(result)
    return result


if __name__ == "__main__":
    main()
