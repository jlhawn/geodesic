#!/usr/bin/env python3
from http import server

"""
Usage:
    python3 httpd.py [port]
"""

class CustomHTTPRequestHandler(server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_custom_headers()
        server.SimpleHTTPRequestHandler.end_headers(self)

    def send_custom_headers(self):
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")

if __name__ == '__main__':
    server.test(HandlerClass=CustomHTTPRequestHandler)
