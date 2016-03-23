from flask import Flask
from app.Vizpymzml.controller import Vizpymzml
from app.Module_Two.controller import WHATEVER

app = Flask(__name__)

app.register_blueprint(Vizpymzml, url_prefix='/')
app.register_blueprint(WHATEVER, url_prefix='/whatever')


