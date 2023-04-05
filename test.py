from tkinter import *
from tkinter import filedialog 
from tkinter import ttk
from configuration import conf
from tkinter import messagebox
from tkinter import font
from multiprocessing.dummy import Process
import threading
import os
import sys
import re
import time

import configuration
import subprocess


p1 = subprocess.Popen(['python','database.py'], shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, bufsize = 1)
