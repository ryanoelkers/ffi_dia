""" This class is used for basic functions not spcecific to this code base such as:
 logging, file writing, and testing."""
from config import Configuration
import logging
import os


class Utils:

    @staticmethod
    def log(statement, level, prnt):
        """ The logger function to log all activity form the program funcationality to both the screen, and file.

        :argument statement: A string which shows want needs to be logged
        :argument level:- The type of statement: info, debug, etc
        :argument prnt: A string with Y or N to either print to the screen and file, or file only

        :return - Nothing is returned, but the log is updated and possibly the log is printed to the screen
        """

        # create the logger
        logger = logging.getLogger()

        if not getattr(logger, 'handler_set', None):
            logger.setLevel(logging.DEBUG)

            # create console handelr and set level to debug
            ch = logging.StreamHandler()
            ch.setLevel(logging.DEBUG)

            # create the formatter
            formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")

            # add formatter to ch
            ch.setFormatter(formatter)

            # add ch to logger
            logger.addHandler(ch)

            # 'set' Handler
            logger.handler_set = True

        if level == 'info':
            logger.info(statement)
        if level == 'debug':
            logger.info(statement)
        if level == 'warning':
            logger.info(statement)
        if level == 'error':
            logger.info(statement)
        if level == 'critical':
            logger.info(statement)

    @staticmethod
    def write_txt(path, typ, line):
        """ This function will either write or append a text file, primarily used for results or logging.

        :parameter path: The location for the file
        :parameter typ: The write type either 'w', or 'a'
        :parameter line: The line you want to write to the file

        :return: Nothing is returned, the file is simply written or appended.

        """

        # open the file for writing
        file = open(path, typ)

        # write the line to the file
        file.write(line)

        # close the file
        file.close()

    @staticmethod
    def get_file_list(path, file_ext):
        """ This function will return the files in a given directory without their extension.
        :parameter path - A string with the path the file.
        :parameter file_ext - The file type to make a list, generally something like *.fits.

        :return files_no_ext - A list of files without the extension."""

        # get the files in the given path directory with the given file extension
        file_list = [f for f in os.listdir(path) if f.endswith(file_ext)]

        # sort based on the number of the image, first taken image will be first
        file_list.sort()

        return file_list

    @staticmethod
    def create_directories(directory_list):
        """ This function will check for each directory in directory list, and create it if it doesn't already exist

        :parameter directory_list - The list of directories to create

        :return - Nothing is returned but directories are created if necessary
        """

        for path in directory_list:
            # if the path does not exist then create it!
            if os.path.exists(path) is False:
                os.mkdir(path)
                Utils.log(path + ' created.', 'info', Configuration.LOG_SCREEN)

    @staticmethod
    def get_index(idx, num, limit):
        """ This function gets the upper and lower index values based on the current index, and the size of a vector

        :parameter idx - The current index
        :parameter num - The number of indeces to find
        :parameter limit - The length of the array

        :return udx, ldx - the upper and lower limits of the current search
        """
        # a basic arithematic to get the limits
        udx = idx + (num / 2)
        ldx = idx - (num / 2)

        # update the upper/lower based on the limits of the vector
        # if the lower limit is below zero, force it to be zero and then bump the upper limit to be the offset
        if (ldx < 0) & (idx < (num / 2)):
            ldx = 0
            udx = num
        # if the upper index is above the limit, then we need to set the upper index as the higher limit
        if udx > limit:
            udx = limit
            ldx = limit - num

        return udx, ldx
