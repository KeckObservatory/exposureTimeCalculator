from flask import Flask, render_template, request, redirect, url_for, session
from etc_gui import etc_gui 

app = Flask(__name__)
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'

# get_resource_as_string is needed to transfer the CSS 
# back to the www/www2 server

def get_resource_as_string(name, charset='utf-8'):
    with app.open_resource(name) as f:
        return f.read().decode(charset)

app.jinja_env.globals['get_resource_as_string'] = get_resource_as_string

# Establish etcgui route

@app.route('/etcgui/',methods=['POST','GET'])
def etcgui():
    return etc_gui()

if __name__ == '__main__':
    host = '0.0.0.0'
    port = 51996
    debug = True
    app.run(host=host,port=port,debug=debug)
