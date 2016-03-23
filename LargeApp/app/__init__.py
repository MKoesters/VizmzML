from flask import Flask, render_template

from app.VizpyML.controllers import Vizpymzml
from app.Module_Two.controllers import WHATEVER

app = Flask(__name__)
app.config.from_object('config')



app.register_blueprint(Vizpymzml, url_prefix='/')
#app.register_blueprint(WHATEVER, url_prefix='/whatever')


