# encoding: utf-8
"""
@author: Xin Zhang
@contact: meetdevin.zh@outlook.com
@file: exception_message.py
@time: 3/6/19
@desc: 自定义异常类
"""


class ExceptionPassing(Exception):
    """
    继承自基类 Exception
    """
    def __init__(self, *message, expression=None):
        super(ExceptionPassing, self).__init__(message)
        self.expression = expression
        self.message = str.join('', [str(a) for a in message])
